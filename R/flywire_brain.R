library(Matrix)
library(dplyr)

new_flywire_brain <- function(mdn, mdn_edges, up_meta = NULL) {

  edges <- mdn_edges |>
    mutate(
      pre_id  = as.character(pre_pt_root_id),
      post_id = as.character(post_pt_root_id)
    )

  upstream_ids <- unique(edges$pre_id)
  mdn_ids <- as.character(mdn$root_id)

  pre_index <- match(edges$pre_id, upstream_ids)
  post_index <- match(edges$post_id, mdn_ids)

  W <- sparseMatrix(
    i = post_index,
    j = pre_index,
    x = edges$syn_count,
    dims = c(
      length(mdn_ids),
      length(upstream_ids)
    ),
    dimnames = list(
      mdn_ids,
      upstream_ids
    )
  )

  list(
    W = W,

    upstream_ids = upstream_ids,
    mdn_ids = mdn_ids,

    upstream_state = setNames(
      numeric(length(upstream_ids)),
      upstream_ids
    ),

    mdn_state = setNames(
      numeric(length(mdn_ids)),
      mdn_ids
    ),

    mdn_side = setNames(
      mdn$side,
      mdn_ids
    ),

    metadata = up_meta
  )
}

flywire_step <- function(brain) {

  x <- brain$W %*% brain$upstream_state

  x <- as.numeric(x)

  names(x) <- brain$mdn_ids

  # Saturating activation.
  #
  # 20 synaptic-equivalent units ~= substantial activation
  # for our simulation.
  x <- 1 - exp(-x / 20)

  brain$mdn_state <- x

  brain
}

mdn_output <- function(brain) {

  left_ids <- names(brain$mdn_side)[
    brain$mdn_side == "left"
  ]

  right_ids <- names(brain$mdn_side)[
    brain$mdn_side == "right"
  ]

  left <- mean(
    brain$mdn_state[left_ids]
  )

  right <- mean(
    brain$mdn_state[right_ids]
  )

  c(
    BACKWARD = mean(c(left, right)),
    LEFT_MDN = left,
    RIGHT_MDN = right,
    LATERALIZATION = left - right
  )
}

build_flywire_network <- function(
    second_edges,
    mdn_edges,
    mdn
) {

  id_character <- function(x) {
    if (inherits(x, "integer64")) bit64::as.character.integer64(x) else as.character(x)
  }

  e2 <- second_edges |>
    transmute(
      pre  = id_character(pre_pt_root_id),
      post = id_character(post_pt_root_id),
      weight = as.numeric(syn_count)
    )

  e1 <- mdn_edges |>
    transmute(
      pre  = if ("pre_id" %in% names(mdn_edges)) pre_id else id_character(pre_pt_root_id),
      post = if ("post_id" %in% names(mdn_edges)) post_id else id_character(post_pt_root_id),
      weight = as.numeric(syn_count)
    )

  edges <- bind_rows(e2, e1) |>
    group_by(pre, post) |>
    summarise(
      weight = sum(weight),
      .groups = "drop"
    )

  neuron_ids <- unique(c(
    edges$pre,
    edges$post
  ))

  index <- setNames(
    seq_along(neuron_ids),
    neuron_ids
  )

  W <- sparseMatrix(
    i = unname(index[edges$post]),
    j = unname(index[edges$pre]),
    x = edges$weight,

    dims = c(
      length(neuron_ids),
      length(neuron_ids)
    ),

    dimnames = list(
      neuron_ids,
      neuron_ids
    )
  )

  mdn_ids <- as.character(mdn$root_id)

  brain <- list(
    W = W,
    edges = edges,

    neuron_ids = neuron_ids,

    state = setNames(
      numeric(length(neuron_ids)),
      neuron_ids
    ),

    mdn_ids = mdn_ids,

    mdn_side = setNames(
      mdn$side,
      mdn_ids
    )
  )

  brain$W_norm <- normalize_connectome(brain$W)
  brain
}

normalize_connectome <- function(W) {

  input_strength <- Matrix::rowSums(W)

  input_strength[
    input_strength == 0
  ] <- 1

  Diagonal(
    x = 1 / input_strength
  ) %*% W
}

flywire_brain_step <- function(
    brain,
    stimulus = NULL,
    decay = 0.2,
    propagation = 0.8
) {

  # Проверяем/восстанавливаем имена state
  if (
    is.null(names(brain$state)) ||
    !identical(names(brain$state), brain$neuron_ids)
  ) {
    names(brain$state) <- brain$neuron_ids
  }

  # Сигнал через реальные связи FlyWire
  recurrent <- as.numeric(
    brain$W_norm %*% brain$state
  )

  names(recurrent) <- brain$neuron_ids

  # Внешний сенсорный вход
  input <- numeric(
    length(brain$state)
  )

  names(input) <- brain$neuron_ids

  if (!is.null(stimulus)) {

    common <- intersect(
      names(stimulus),
      brain$neuron_ids
    )

    input[common] <- stimulus[common]
  }

  # Динамика
  x <-
    decay * brain$state +
    propagation * recurrent +
    input

  # Ограничиваем 0...1
  new_state <- pmin(
    1,
    pmax(0, x)
  )

  # КРИТИЧЕСКИ ВАЖНО:
  # возвращаем FlyWire ID в качестве names
  names(new_state) <- brain$neuron_ids

  brain$state <- new_state

  brain
}

make_visual_stimulus <- function(
    brain,
    intensity = 1,
    group = "LT51"
) {
  ids <- intersect(brain$visual_groups[[group]], brain$neuron_ids)

  if (!length(ids)) {
    return(setNames(numeric(), character()))
  }

  setNames(
    rep(intensity, length(ids)),
    ids
  )
}


run_visual_test <- function(
    brain,
    group,
    steps = 15,
    intensity = 1
) {
  brain$state[] <- 0

  result <- vector("list", steps)

  for (t in seq_len(steps)) {

    stimulus <- if (t == 1) {
      make_visual_stimulus(
        brain,
        intensity = intensity,
        group = group
      )
    } else {
      NULL
    }

    brain <- flywire_brain_step(
      brain,
      stimulus = stimulus
    )

    result[[t]] <- c(
      time = t,
      read_mdn(brain)
    )
  }

  x <- as.data.frame(
    do.call(rbind, result)
  )

  x$group <- group

  x
}

read_mdn <- function(brain) {

  mdn_ids <- as.character(
    brain$mdn_ids
  )

  if (is.null(names(brain$state))) {
    stop("brain$state не имеет names")
  }

  missing <- setdiff(
    mdn_ids,
    names(brain$state)
  )

  if (length(missing) > 0) {
    stop(
      "MDN отсутствуют в brain$state: ",
      paste(missing, collapse = ", ")
    )
  }

  activity <- brain$state[mdn_ids]

  sides <- as.character(
    brain$mdn_side[mdn_ids]
  )

  left <- activity[
    sides == "left"
  ]

  right <- activity[
    sides == "right"
  ]

  safe_mean <- function(x) if (length(x)) mean(x) else 0

  c(
    BACKWARD = safe_mean(activity),
    LEFT_MDN = safe_mean(left),
    RIGHT_MDN = safe_mean(right),
    LATERALIZATION =
      safe_mean(left) - safe_mean(right)
  )
}
