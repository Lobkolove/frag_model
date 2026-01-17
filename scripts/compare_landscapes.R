# Script to compare the output of the old function used for landscape generation
# so far (nlm_fbm) and the new proposed one (fbm_fft), which has the benefit of
# generating toroidal grids.

# Note that output from the old function needed to be generated in a separate
# session running an older R version due to dependency issues.
library(viridis)
source("Model/src/landscape.R")

# Import landscape from old function:
old_ls <- raster::raster("data-raw/nlm_fbm_landscape.tif")

# Generate a landscape with new function using the same parameters:
new_ls <- fbm_fft(gr_size = 128, ac_amount = .7, seed = 42)
new_ls


# 0. Quick structural sanity check ----------------------------------------
raster::compareRaster(old_ls, new_ls)


# 1. Visual Comparison ----------------------------------------------------

par(mfrow = c(1, 2))
image(old_ls, asp = 1, main = "nlm_fbm landscape", col = viridis(100))
image(new_ls, asp = 1, main = "fbm_fft landscape", col = viridis(100))
par(mfrow = c(1, 1))


# 2. Compare distributions of values --------------------------------------

par(mfrow = c(1, 2))
hist(values(old_ls), breaks = 50, probability = TRUE, 
     col = scales::alpha("midnightblue", 0.5), 
     main = "Old (RandomFields)", xlab = "Value")
hist(values(new_ls), breaks = 50, probability = TRUE, 
     col = scales::alpha("violetred4", 0.5), 
     main = "New (FFT)", xlab = "Value")
par(mfrow = c(1, 1))

# Overlay for direct comparison
hist(values(old_ls), breaks = 50, probability = TRUE, 
     col = scales::alpha("midnightblue", 0.5), 
     ylim = c(0, max(density(values(old_ls))$y, 
                     density(values(new_ls))$y) * 1.1), 
     main = "Distribution comparison", xlab = "Value")
hist(values(new_ls), breaks = 50, probability = TRUE, 
     col = scales::alpha("violetred4", 0.3), add = TRUE)

