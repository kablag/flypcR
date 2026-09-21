apply_action <- function(world, action) {

  p <- world$pipette

  if (action == "LEFT")
    p$x <- max(1, p$x - 1)

  if (action == "RIGHT")
    p$x <- min(12, p$x + 1)

  if (action == "UP")
    p$y <- min(8, p$y + 1)

  if (action == "DOWN")
    p$y <- max(1, p$y - 1)

  if (action == "DISPENSE") {

    i <- which(
      world$wells$col == p$x &
      match(world$wells$row, LETTERS) == p$y
    )

    # A completed well stays completed on subsequent UI clicks.
    world$wells$volume[i] <- pmin(
      p$volume,
      world$wells$volume[i] + p$volume
    )
  }

  world$pipette <- p

  world
}

reward <- function(world) {

  target <- world$wells[
    world$wells$target,
  ]

  if (target$volume == 10)
    return(100)

  target_xy <- well_xy(
    target$row,
    target$col
  )

  distance <-
    abs(world$pipette$x - target_xy["x"]) +
    abs(world$pipette$y - target_xy["y"])

  -distance
}

