# Demonstration of use of toroidal_clump function to identify patches in a raster with toroidal (wrap-around) boundaries.
library(raster)
library(viridis)
library(Polychrome)
source("Model/src/landscape.R")
source("Model/src/fragmentation.R")
source("R/toroidal_clump.R")

seed_landscape <- 42
seed_fragment <- 43

# Create landscape
full_ls <- fbm_fft(gr_size = 128, seed = seed_landscape)
par(mar = c(0.2, 0.2, 0.2, 0.2))
image(full_ls, col = viridis(100), asp = 1, axes = FALSE, xlab = "", ylab = "")

# Fragment landscape
fragmented_ls <- ls_mask(full_ls, habitat = 0.15, fragmentation = 0.4, seed = seed_fragment)
image(fragmented_ls, col = viridis(100), asp = 1, axes = FALSE, xlab = "", ylab = "")

# Normal clumping (no toroidal)
rclumped <- raster::clump(fragmented_ls, directions = 4, gaps = FALSE)
r_ids <- unique(na.omit(raster::getValues(rclumped)))
image(rclumped, col = Polychrome::palette36.colors(length(r_ids)), 
      asp = 1, xlab = "", ylab = "")

# Identify patches with toroidal clumping
tclumped <- toroidal_clump(fragmented_ls, directions = 4)
t_ids <- na.omit(unique(raster::getValues(tclumped)))
pal <- Polychrome::palette36.colors(length(t_ids))
image(tclumped, col = pal, 
      asp = 1, xlab = "", ylab = "")
for (id in t_ids) {
  tmp <- tclumped
  tmp[tmp != id]  <- NA
  image(tmp, main = paste("Patch ID:", id), asp = 1, col = pal[id])
}


# Test with plain matrix
alt_ls <- fbm_fft(gr_size = 128, seed = seed_landscape, raster = FALSE)
alt_fragmented <- ls_mask(alt_ls, habitat = 0.15, fragmentation = 0.4, seed = seed_fragment)
t2 <- toroidal_clump(alt_fragmented, directions = 4)
image(t2, col = pal, axes = F)
