# Script to test outputs of the functions simulating the fragmentation step.
library(raster)
library(data.table)
library(viridis)
source("Model/parameters.R")
source("Model/src/initialize.R")
source("Model/src/generate_grid.R")
source("Model/src/landscape.R")
source("Model/src/generate_agents.R")
source("Model/src/distribute_agents.R")

source("Model/src/cookie_cutting.R")
source("Model/src/ls_mask.R")

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


# Initialize list for fragmented grids
frag_seq <- list(cookie_cutter = list(),
                 mask = list())

# Generate fragmented grids
for (i in 1:length(frag_levels)) {
  
  # With cookie cutter function
  full_ck <- cookie_cutting(
    grid = full_grid,
    agents = agents,
    agents_grid = agents_grid,
    habitat = habitat_ratio,
    fragmentation = frag_levels[i]
  )
  frag_seq$cookie_cutter[i] <- full_ck$grid
  
  # With mask function
  full_ms <- ls_mask(
    grid = full_grid,
    agents = agents,
    agents_grid = agents_grid,
    habitat = habitat_ratio,
    fragmentation = frag_levels[i]
  )
  frag_seq$mask[i] <- full_ms$grid
}

# Visualise fragmented grids
par(mfrow = c(1,2), pty = "s")
for (i in 1:length(frag_levels)) {
  
  par(mar = c(2,2,2,1))
  image(frag_seq$cookie_cutter[[i]], main = "Cookie Cutter",
        asp = 1, col = viridis(100),
        xaxt = "n", yaxt = "n")
  par(mar = c(2,1,2,2))
  image(frag_seq$mask[[i]], main = "Mask",
        asp = 1, col = viridis(100),
        xaxt = "n", yaxt = "n")
  mtext(paste0("Fragmentation = ", frag_levels[i]), 
        side = 3, line = -3, cex = 1.2, outer = TRUE)
}


