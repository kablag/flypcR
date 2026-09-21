load_flywire_data <- function(
    file = "data/flywire_subconnectome_v783.rds"
) {
  if (!file.exists(file)) {
    stop("FlyWire data file not found: ", file)
  }

  data <- readRDS(file)
  required <- c("mdn", "mdn_edges", "second_edges", "visual_groups")
  missing <- setdiff(required, names(data))

  if (length(missing)) {
    stop("Invalid FlyWire data; missing: ", paste(missing, collapse = ", "))
  }

  data
}
