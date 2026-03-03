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
par(mar = c(0, 0, 0, 0))
image(full_ls, col = viridis(100), asp = 1, axes = FALSE, xlab = "", ylab = "")

# Fragment landscape
fragmented_ls <- ls_mask(full_ls, habitat = 0.15, fragmentation = 0.3, seed = seed_fragment)
image(fragmented_ls, col = viridis(100), asp = 1, axes = FALSE, xlab = "", ylab = "")
image(is.na(fragmented_ls), col = c(NA, "grey90"), asp = 1, add = TRUE)

# Export fragmented landscape
png("pics/fragmented_ls.png", width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(fragmented_ls, col = viridis(100), asp = 1, axes = FALSE)
image(is.na(fragmented_ls), col = c(NA, "grey90"), asp = 1, add = TRUE)
dev.off()

# Normal clumping (no toroidal)
rclumped <- raster::clump(fragmented_ls, directions = 4, gaps = FALSE)
r_ids <- unique(na.omit(raster::getValues(rclumped)))
pal <- as.character(paletteer::paletteer_d("MoMAColors::Sidhu"))
image(rclumped, col = pal[1:length(r_ids)], 
      asp = 1, axes = FALSE)
image(is.na(fragmented_ls), col = c(NA, "grey90"), asp = 1, add = TRUE)
box(col = "black", lwd = 1)

# Export to file
png("pics/r_clumps.png", width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(rclumped, col = pal[1:length(r_ids)], 
      asp = 1, axes = FALSE, xlab = "", ylab = "")
image(is.na(rclumped), col = c(NA, "grey90"), asp = 1, add = TRUE)
box(col = "black", lwd = 1)
dev.off()

# Identify patches with toroidal clumping
tclumped <- toroidal_clump(fragmented_ls, directions = 4)
t_ids <- na.omit(unique(raster::getValues(tclumped)))
par(mar = c(0, 0, 0, 0))
image(tclumped, col = pal[1:length(t_ids)], 
      asp = 1, axes = FALSE)
image(is.na(tclumped), col = c(NA, "grey90"), asp = 1, add = TRUE)
box(col = "black", lwd = 1)

# Export to file
png("pics/t_clumps.png", width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(tclumped, col = pal[1:length(t_ids)], 
      asp = 1, axes = FALSE, xlab = "", ylab = "")
image(is.na(tclumped), col = c(NA, "grey90"), asp = 1, add = TRUE)
box(col = "black", lwd = 1)
dev.off()


# Test with plain matrix
alt_ls <- fbm_fft(gr_size = 128, seed = seed_landscape, raster = FALSE)
alt_fragmented <- ls_mask(alt_ls, habitat = 0.15, fragmentation = 0.4, seed = seed_fragment)
t2 <- toroidal_clump(alt_fragmented, directions = 4)
image(t2, col = pal, axes = F)
