new_world <- function() {

  wells <- expand.grid(
    row = LETTERS[1:8],
    col = 1:12
  )

  wells$volume <- 0
  wells$target <- FALSE

  # Первая задача — попасть в A1
  wells$target[wells$row == "A" & wells$col == 1] <- TRUE

  list(
    pipette = list(
      x = 6,
      y = 4,
      volume = 10
    ),
    wells = wells
  )
}


well_xy <- function(row, col) {
  c(
    x = col,
    y = match(row, LETTERS[1:8])
  )
}


observe_world <- function(world) {

  target <- world$wells[world$wells$target, ]

  target_xy <- well_xy(
    target$row,
    target$col
  )

  dx <- target_xy["x"] - world$pipette$x
  dy <- target_xy["y"] - world$pipette$y

  c(
    LEFT  = as.numeric(dx < 0),
    RIGHT = as.numeric(dx > 0),
    UP    = as.numeric(dy > 0),
    DOWN  = as.numeric(dy < 0),
    HERE  = as.numeric(dx == 0 && dy == 0)
  )
}