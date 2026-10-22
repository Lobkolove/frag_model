library(dplyr)
library(ggplot2)
library(testthat)
source("R/sample_cells.R")
source("R/toroidal_clump.R")

source("R/toroidal_dist.R")

# Example coords from saved state
states <- readRDS("data-raw/states_frag_0.75_sim_6.rds")
state1 <- states$post_fragmentation

samp <- sample_cells(state1, method = "random", n_samples = 30)

coords <- samp |> 
  select(x_loc, y_loc)

# dist object
D <- stats::dist(coords)
D
length(D)

grid_size <- samp$grid_size[1]
T <- toroidal_dist(coords, grid_size)
T
length(T)
class(T)
names(T)

test_that("toroidal_dist matches stats::dist on small grid", {
  # Small grid where no wrapping occurs
  coords <- matrix(c(0, 0, 5, 5, 10, 10), ncol = 2, byrow = TRUE)
  grid_size <- 100
  
  torus_d <- toroidal_dist(coords, grid_size)
  euclidean_d <- stats::dist(coords)
  
  expect_equal(as.numeric(torus_d), as.numeric(euclidean_d), tolerance = 1e-10)
})

test_that("toroidal wrapping is symmetric", {
  # Two points near opposite edges
  coords <- matrix(c(1, 1, 99, 99), ncol = 2, byrow = TRUE)
  grid_size <- 100
  
  d <- as.numeric(toroidal_dist(coords, grid_size))
  
  # Distance should be sqrt(2)*2, not sqrt(2)*98
  expect_lt(d, sqrt(2) * 10)  # Should wrap
})


