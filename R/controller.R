pipette_catalog <- data.frame(
  pipette = c("P1000", "P200", "P20", "P10"),
  min = c(100, 20, 2, 1),
  max = c(1000, 200, 20, 10)
)

choose_pipette <- function(volume) {
  candidates <- lapply(seq_len(nrow(pipette_catalog)), function(i) {
    p <- pipette_catalog[i, ]
    strokes <- ceiling(volume / p$max)
    aliquot <- volume / strokes
    if (aliquot < p$min) return(NULL)
    data.frame(pipette = p$pipette, strokes = strokes, aliquot = aliquot,
               capacity = p$max)
  })
  candidates <- do.call(rbind, candidates)
  if (is.null(candidates) || !nrow(candidates)) stop("Нельзя перенести объём: ", volume)
  candidates <- candidates[order(candidates$strokes, candidates$capacity), ]
  candidates[1, ]
}

build_protocol <- function(n_samples) {
  wells <- pcr_layout(n_samples)
  n <- nrow(wells)
  components <- data.frame(
    stage = "Приготовление мастер-микса",
    source = c("PCR mix 10X", "Праймеры 5X", "Вода"),
    destination = "Мастер-микс 1,5 мл",
    volume = n * c(2.5, 5, 15.5),
    row = NA_character_, col = NA_integer_, new_tip = TRUE
  )
  dispensing <- data.frame(
    stage = "Раздача мастер-микса", source = "Мастер-микс 1,5 мл",
    destination = paste0(wells$row, wells$col), volume = 23,
    row = wells$row, col = wells$col, new_tip = c(TRUE, rep(FALSE, n - 1L))
  )
  templates <- data.frame(
    stage = "Добавление матрицы", source = wells$sample,
    destination = paste0(wells$row, wells$col), volume = 2,
    row = wells$row, col = wells$col, new_tip = TRUE
  )
  tasks <- rbind(components, dispensing, templates)
  choice <- lapply(tasks$volume, choose_pipette)
  tasks$pipette <- vapply(choice, function(x) x$pipette, character(1))
  tasks$strokes <- vapply(choice, function(x) x$strokes, numeric(1))
  tasks$aliquot <- vapply(choice, function(x) x$aliquot, numeric(1))
  tasks
}

deck_position <- function(name, wells) {
  fixed <- data.frame(
    name = c("PCR mix 10X", "Праймеры 5X", "Вода", "Мастер-микс 1,5 мл",
             "P1000", "P200", "P20", "P10",
             "tips_P1000", "tips_P200", "tips_P20", "tips_P10", "Отходы"),
    x = c(1, 3, 5, 7, 1, 3, 5, 7, 10, 12, 14, 16, 19),
    y = c(2, 2, 2, 2, rep(15, 9)), stringsAsFactors = FALSE
  )
  hit <- fixed[fixed$name == name, ]
  if (nrow(hit)) return(c(x = hit$x[1], y = hit$y[1]))
  well <- wells[paste0(wells$row, wells$col) == name, ]
  if (nrow(well)) return(c(x = 10 + well$col[1] - 1, y = 2 + match(well$row[1], LETTERS) - 1))
  sample <- match(name, wells$sample)
  if (!is.na(sample)) return(c(x = 1 + ((sample - 1) %% 7), y = 5 + ((sample - 1) %/% 7)))
  stop("Неизвестный объект на столе: ", name)
}

build_run_actions <- function(tasks, wells) {
  out <- list(); add <- function(target, verb, task, volume = 0) {
    xy <- deck_position(target, wells)
    out[[length(out) + 1L]] <<- data.frame(
      stage = task$stage, verb = verb, target = target, x = xy[["x"]], y = xy[["y"]],
      pipette = task$pipette, volume = volume, task = as.integer(task$.id),
      stringsAsFactors = FALSE)
  }
  tasks$.id <- seq_len(nrow(tasks)); current_pipette <- NA_character_; has_tip <- FALSE
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, ]
    if (!identical(current_pipette, task$pipette)) {
      if (has_tip) { add("Отходы", "Сбросить наконечник", task); has_tip <- FALSE }
      add(task$pipette, "Взять пипетку", task); current_pipette <- task$pipette
    }
    if (task$new_tip) {
      if (has_tip) add("Отходы", "Сбросить наконечник", task)
      add(paste0("tips_", task$pipette), "Надеть наконечник", task); has_tip <- TRUE
    }
    for (stroke in seq_len(task$strokes)) {
      add(task$source, "Набрать", task, task$aliquot)
      add(task$destination, "Выдать", task, task$aliquot)
    }
  }
  if (has_tip) add("Отходы", "Сбросить наконечник", tasks[nrow(tasks), ])
  do.call(rbind, out)
}

apply_action <- function(world, action = NULL) {
  if (world$action_index >= nrow(world$actions)) return(world)
  world$step_count <- world$step_count + 1L
  a <- world$actions[world$action_index + 1L, ]
  dx <- a$x - world$fly$x; dy <- a$y - world$fly$y
  if (dx != 0 || dy != 0) {
    if (abs(dx) >= abs(dy) && dx != 0) world$fly$x <- world$fly$x + sign(dx)
    else world$fly$y <- world$fly$y + sign(dy)
    return(world)
  }

  if (a$verb == "Взять пипетку") world$mounted_pipette <- a$pipette
  if (a$verb == "Надеть наконечник") {
    world$has_tip <- TRUE; world$tips[a$pipette] <- world$tips[a$pipette] + 1L
  }
  if (a$verb == "Сбросить наконечник") world$has_tip <- FALSE
  if (a$verb == "Набрать") world$carried_volume <- a$volume
  if (a$verb == "Выдать") {
    task <- world$tasks[a$task, ]
    if (task$stage == "Приготовление мастер-микса") {
      world$master[task$source] <- world$master[task$source] + a$volume
      world$master_volume <- world$master_volume + a$volume
    } else {
      well <- which(world$wells$row == task$row & world$wells$col == task$col)
      if (task$stage == "Раздача мастер-микса") {
        world$wells$master_mix[well] <- world$wells$master_mix[well] + a$volume
        world$master_volume <- world$master_volume - a$volume
      } else world$wells$template[well] <- world$wells$template[well] + a$volume
      world$wells$volume[well] <- world$wells$master_mix[well] + world$wells$template[well]
    }
    world$carried_volume <- 0
  }
  world$log <- c(world$log, sprintf("%s: %s %s%s", a$stage, a$verb, a$target,
    if (a$volume > 0) paste0(" — ", round(a$volume, 2), " мкл") else ""))
  world$action_index <- world$action_index + 1L
  world$task_index <- if (world$action_index) max(world$actions$task[seq_len(world$action_index)]) else 0L
  world
}
