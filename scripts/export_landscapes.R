library(raster)
library(viridis)
library(ggplot2)
source("Model/src/landscape.R")
source("Model/src/fragmentation.R")
source("R/tile_raster.R")

# Set parameters
gr_size <- 400
ac_amount <- 0.5

# Unique identifier (change for new images)
id <- 3

# Full landscape ---------------------------------------------------------

seed_landscape <- 42

fbm <- fbm_fft(gr_size = gr_size, ac_amount = ac_amount, raster = T, seed = seed_landscape)
# As raster to have consistent plotting with ls_mask output

# Plot to screen
par(mar = c(0, 0, 0, 0))
image(fbm, asp = 1, axes = FALSE,
      col = viridis(100))
dev.off()

# Export to file
filename <- paste0("pics/fbm_full_ac_", ac_amount, "_", id, ".png")
png(filename, width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(fbm, asp = 1, axes = FALSE,
      col = viridis(100))
dev.off()


# Tiled landscape --------------------------------------------------------

tiled <- tile_raster(fbm, n = 3)

# Plot to screen
par(mar = c(0, 0, 0, 0))
image(tiled, asp = 1, axes = FALSE,
      col = viridis(100))
dev.off()

# Export to file
filename <- paste0("pics/fbm_tiled_ac_", ac_amount, "_", id, ".png")
png(filename, width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(tiled, asp = 1, axes = FALSE,
      col = viridis(100))
dev.off()


# Fragmented landscape ---------------------------------------------------

fragmentation <- 0.5
seed_fragment <- seed_landscape + 1 
# Just to have a different seed for fragmentation

frag <- ls_mask(fbm, habitat = 0.15, fragmentation = fragmentation, seed = seed_fragment)

# Plot to screen
par(mar = c(0, 0, 0, 0))
image(frag, asp = 1, axes = FALSE,
      col = viridis(100))
# Add masked areas in grey
image(is.na(frag), asp = 1, axes = FALSE,
      col = c(NA, "grey90"), add = TRUE)
dev.off()

# Export to file
filename <- paste0("pics/fbm_fragmented_ac_", ac_amount, "_frag_", fragmentation, "_", id, ".png")
png(filename, width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(frag, asp = 1, axes = FALSE,
      col = viridis(100))
# Add masked areas in grey
image(is.na(frag), asp = 1, axes = FALSE,
      col = c(NA, "grey90"), add = TRUE)
dev.off()
