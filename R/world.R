pcr_layout <- function(n_samples) {
  n_reactions <- as.integer(n_samples) + 5L
  if (n_samples < 0 || n_reactions > 65L) {
    stop("Число реакций должно быть от 5 до 65 для одной пробирки 1,5 мл")
  }

  labels <- c(
    if (n_samples) paste0("ДНК ", seq_len(n_samples)),
    "ПКО WT", "ПКО WT/M", "ПКО M", "ОКО-В", "ОКО (вода)"
  )
  positions <- expand.grid(row = LETTERS[1:8], col = 1:12)
  positions <- positions[order(match(positions$row, LETTERS), positions$col), ]
  positions <- positions[seq_len(n_reactions), ]

  data.frame(
    reaction = seq_len(n_reactions),
    row = as.character(positions$row),
    col = positions$col,
    sample = labels,
    volume = 0,
    master_mix = 0,
    template = 0,
    target = TRUE,
    stringsAsFactors = FALSE
  )
}

new_world <- function(n_samples = 1L) {
  wells <- pcr_layout(n_samples)
  tasks <- build_protocol(n_samples)
  list(
    n_samples = as.integer(n_samples),
    n_reactions = nrow(wells),
    fly = list(x = 1, y = 15),
    mounted_pipette = NA_character_,
    has_tip = FALSE,
    carried_volume = 0,
    wells = wells,
    master = c(`PCR mix 10X` = 0, `Праймеры 5X` = 0, Вода = 0),
    master_volume = 0,
    tasks = tasks,
    task_index = 0L,
    actions = build_run_actions(tasks, wells),
    action_index = 0L,
    step_count = 0L,
    tips = c(P1000 = 0L, P200 = 0L, P20 = 0L, P10 = 0L),
    log = character()
  )
}

well_xy <- function(row, col) c(x = col, y = match(row, LETTERS[1:8]))

observe_world <- function(world) {
  if (world$action_index >= nrow(world$actions)) {
    return(c(LEFT = 0, RIGHT = 0, UP = 0, DOWN = 0, HERE = 1))
  }
  action <- world$actions[world$action_index + 1L, ]
  dx <- action$x - world$fly$x
  dy <- action$y - world$fly$y
  c(LEFT = as.numeric(dx < 0), RIGHT = as.numeric(dx > 0),
    UP = as.numeric(dy > 0), DOWN = as.numeric(dy < 0),
    HERE = as.numeric(dx == 0 && dy == 0))
}

reward <- function(world) {
  if (!nrow(world$wells)) return(0)
  sum(abs(world$wells$master_mix - 23) < 1e-8 &
        abs(world$wells$template - 2) < 1e-8)
}
