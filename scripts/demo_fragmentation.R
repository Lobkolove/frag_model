# Script to test outputs of the functions simulating the fragmentation step.
library(raster)
library(data.table)
library(dplyr)
library(purrr)
library(viridis)
source("Model/parameters.R")
source("Model/src/initialize.R")
source("Model/src/generate_grid.R")
source("Model/src/landscape.R")
source("Model/src/generate_agents.R")
source("Model/src/distribute_agents.R")

source("Model/src/cookie_cutting.R")
source("Model/src/fragmentation.R")


# 0. Setup ----------------------------------------------------------------

# Example parameters:
spatial_ac <- 0.7
niche_breadth <- 0.1
habitat_ratio <- 0.15

# Fragmentation levels
frag_levels <- seq(0.1, 0.9, by = 0.1)

# Initialize model with 100% habitat
model_start <- initialize(
  frag = 0,
  hab  = 1,
  ac   = spatial_ac,
  nb   = niche_breadth
)

# extract simulation space, agents grid and agents list
full_grid    <- model_start$grid
agents_grid  <- model_start$agents_grid
agents       <- model_start$agents

# Visualise full habitat grid
raster::image(full_grid, asp = 1, col = viridis(100))


# 1. Generate outputs -----------------------------------------------------

# Initialise results data frame (one row per frag_level)
results <- tibble(frag_level = frag_levels,
                  cookie_cutter = list(NULL),
                  mask = list(NULL))

# Sequential generation
for (i in seq_along(frag_levels)) {

  # Store all outputs (grid, agents, agents_grid)
  
  results$cookie_cutter[[i]] <- cookie_cutting(grid = full_grid, 
                                               agents = agents, 
                                               agents_grid = agents_grid,
                                               habitat = habitat_ratio, 
                                               fragmentation = frag_levels[i])
  
  results$mask[[i]] <- fragment(grid = full_grid, 
                                agents = agents, 
                                agents_grid = agents_grid,
                                habitat = habitat_ratio, 
                                fragmentation = frag_levels[i],
                                seed = seed)
}


# 2. Compare outputs ------------------------------------------------------

  ## 2.1 Compare grids ####

par(mfrow = c(1,2), pty = "s")
for (i in seq_along(frag_levels)) {
  
  cat("Fragmentation level: ", frag_levels[i], "\n")
  
  old_grid <- results$cookie_cutter[[i]]$grid
  new_grid <- results$mask[[i]]$grid
  
  # Grid geometry check
  if (!raster::compareRaster(old_grid, new_grid, res = TRUE)) {
    cat("Raster comparison failed!", frag_levels[i], "\n")
  }
  
  # Grid values check  
  diff_grid <- old_grid != new_grid
  cat("Grid cells differing:", sum(getValues(diff_grid), na.rm=TRUE), "\n\n")
  
  par(mar = c(2,2,2,1))
  image(old_grid, main = "Cookie Cutter",
        asp = 1, col = viridis(100),
        xaxt = "n", yaxt = "n")
  par(mar = c(2,1,2,2))
  image(new_grid, main = "Mask",
        asp = 1, col = viridis(100),
        xaxt = "n", yaxt = "n")
  mtext(paste0("Fragmentation = ", frag_levels[i]), 
        side = 3, line = -3, cex = 1.2, outer = TRUE)
}

  ## 2.2 Compare agents ####

# ...

  ## 2.3 Compare agents_grid ####

# ...


