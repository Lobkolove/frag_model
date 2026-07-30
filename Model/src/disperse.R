######################## Dispersal function ########################################

# The function takes the current location of an individual and returns a new location
# based on a log-normal dispersal kernel or an exponential dispersal kernel

disperse <- function(cur_loc, d_sd, d_mean) {
  dis_sd <- d_sd
  dis_mean <- d_mean

  if (switch$kernel_type == 0) {
    sigma <- sqrt(log(1 + (dis_sd * dis_sd) / (dis_mean * dis_mean)))
    mu <- log(dis_mean) - 0.5 * sigma * sigma

    distance <- rlnorm(1, mu, sigma)
  } else if (switch$kernel_type == 1) {
    distance <- rexp(1, 1 / dis_mean)
  } else {
    print("please check switches")
  }

  direction <- runif(1, min = 0, max = 2 * pi)

  new_x <- cur_loc[1] + round(cos(direction) * distance)
  new_y <- cur_loc[2] + round(sin(direction) * distance)

  new_loc <- c(new_x, new_y)
  return(new_loc)
}


######################## Toroidal dispersal ########################################

# Wraps coordinates on a toroidal grid: moves that exit one side re-enter from
# the opposite side. Keeps the original kernel logic but enforces modular
# coordinates instead of rejecting off-grid moves.
toroidal_disperse <- function(
  cur_loc,
  d_sd,
  d_mean,
  grid_size,
  kernel = c("exponential", "log-normal")
) {
  start.time <- Sys.time()

  kernel <- match.arg(kernel)

  dis_sd <- d_sd
  dis_mean <- d_mean

  if (kernel == "log-normal") {
    sigma <- sqrt(log(1 + (dis_sd * dis_sd) / (dis_mean * dis_mean)))
    mu <- log(dis_mean) - 0.5 * sigma * sigma
    distance <- rlnorm(1, mu, sigma)
  } else if (kernel == "exponential") {
    distance <- rexp(1, 1 / dis_mean)
  } else {
    print("please check switches")
  }

  direction <- runif(1, min = 0, max = 2 * pi)

  raw_x <- cur_loc[1] + round(cos(direction) * distance)
  raw_y <- cur_loc[2] + round(sin(direction) * distance)

  # wrap to 1..grid_size
  new_x <- ((raw_x - 1) %% grid_size) + 1
  new_y <- ((raw_y - 1) %% grid_size) + 1

  new_loc <- c(new_x, new_y)

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  # cat("Time taken: ", time.taken)

  return(new_loc)
}

random_disperse <- function(
  nrow = 50,
  ncol = nrow,
  grid = NULL,
  force_habitat = FALSE,
  habitat_cells = NULL,
  seed = NULL
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (force_habitat && is.null(grid) && is.null(habitat_cells)) {
    stop("Either a grid or habitat_cells must be provided when force_habitat is TRUE.")
  }

  if (!is.null(grid)) {
    nrow <- nrow(grid)
    ncol <- ncol(grid)
  }

  if (force_habitat) {
    if (is.null(habitat_cells)) {
      habitat_cells <- which(!is.na(raster::getValues(grid)))
    }
    destination <- sample(habitat_cells, 1)
    new_loc <- c(ceiling(destination / ncol), (destination - 1) %% ncol + 1)
  } else {
    new_loc <- c(sample(nrow, 1), sample(ncol, 1))
  }

  return(new_loc)
}
