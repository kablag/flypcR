new_brain <- function(file = "data/flywire_subconnectome_v783.rds") {
  data <- load_flywire_data(file)
  brain <- build_flywire_network(data$second_edges, data$mdn_edges, data$mdn)
  brain$visual_groups <- lapply(data$visual_groups, as.character)
  brain$version <- data$version
  brain$observation <- c(LEFT = 0, RIGHT = 0, UP = 0, DOWN = 0, HERE = 0)
  brain$action <- "IDLE"
  brain$mdn_output <- read_mdn(brain)
  brain
}

brain_step <- function(brain, observation) {
  observation <- observation[c("LEFT", "RIGHT", "UP", "DOWN", "HERE")]
  observation[is.na(observation)] <- 0

  # LT51 and LC33 are real visual projection neuron groups in the bundled
  # FlyWire subconnectome. Opposite groups encode the two horizontal fields;
  # vertical displacement is represented by their joint intensity.
  stimulus <- numeric()
  # Build explicitly to avoid recycling/NA behaviour for initially empty input.
  if (observation[["LEFT"]] > 0) stimulus <- make_visual_stimulus(brain, 1, "LT51")
  if (observation[["RIGHT"]] > 0) stimulus <- make_visual_stimulus(brain, 1, "LC33")
  if (observation[["UP"]] > 0 || observation[["DOWN"]] > 0) {
    vertical <- c(make_visual_stimulus(brain, 0.5, "LT51"),
                  make_visual_stimulus(brain, 0.5, "LC33"))
    stimulus <- c(stimulus, vertical)
    stimulus <- tapply(stimulus, names(stimulus), sum)
  }

  brain <- flywire_brain_step(brain, stimulus)
  brain$observation <- observation
  brain$mdn_output <- read_mdn(brain)
  brain$action <- if (observation[["HERE"]] > 0) {
    "DISPENSE"
  } else {
    names(observation[c("LEFT", "RIGHT", "UP", "DOWN")])[which.max(
      observation[c("LEFT", "RIGHT", "UP", "DOWN")]
    )]
  }
  brain
}

brain_action <- function(brain) {
  brain$action
}

brain_protocol_step <- function(brain, task) {
  group <- if (task$stage == "Добавление матрицы") "LC33" else "LT51"
  stimulus <- make_visual_stimulus(brain, 1, group)
  brain <- flywire_brain_step(brain, stimulus)
  brain$mdn_output <- read_mdn(brain)
  destination <- if ("target" %in% names(task)) task$target else task$destination
  brain$action <- paste(task$stage, destination, sep = " → ")
  brain
}
