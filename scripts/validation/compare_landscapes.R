# Script to compare the output of the old function used for landscape generation
# so far (nlm_fbm) and the new proposed one (fbm_fft), which has the benefit of
# generating toroidal grids.

# Note that output from the old function needed to be generated in a separate
# session running an older R version due to dependency issues.
library(tidyr)
library(ggplot2)
library(viridis)
source("Model/src/landscape.R")

# Import landscapes from old function into a raster stack:
old_files <- list.files(path = "data-raw/landscapes",
                        pattern = "nlm_fbm_ac.*\\.tif$",
                        full.names = TRUE)
old_stack <- raster::stack(old_files)
names(old_stack)

# Generate landscapes with new function using the same parameters:
gr_size <- 128
ac_vals <- seq(0.01, 0.91, 0.10)
resolution <- 1
set.seed(42)

new_stack <- raster::stack()
for (ac_amount in ac_vals) {
  new_ls <- fbm_fft(gr_size, resolution, ac_amount)
  names(new_ls) <- paste0("fbm_fft_ac", ac_amount)
  new_stack <- raster::addLayer(new_stack, new_ls)
  rm(new_ls)
}

raster::nlayers(new_stack)
names(new_stack)
raster::plot(new_stack$fbm_fft_ac0.71, col = viridis(100))


# 0. Quick structural sanity check ----------------------------------------

for (ac_amount in ac_vals) {
  old_ls <- old_stack[[match(ac_amount, ac_vals)]]
  new_ls <- new_stack[[match(ac_amount, ac_vals)]]
  
  if (raster::compareRaster(old_ls, new_ls)) 
    cat("Sanity check passed for ac_amount = ", ac_amount, "\n")
}


# 1. Visual Comparison ----------------------------------------------------

par(mfrow = c(1, 2))
for (ac_amount in ac_vals) {
  old_ls <- old_stack[[match(ac_amount, ac_vals)]]
  new_ls <- new_stack[[match(ac_amount, ac_vals)]]
  
  raster::image(old_ls, asp = 1, axes = F, xlab = "",
                sub = "nlm_fbm", col = viridis(100))
  raster::image(new_ls, asp = 1, axes = F, xlab = "",
                sub = "fbm_fft", col = viridis(100))
  mtext(paste0("ac_amount = ", ac_amount), 
        side = 3, line = -2, outer = TRUE)
}
par(mfrow = c(1, 1))


# 2. Compare distributions of values --------------------------------------

# Histograms
for (ac_amount in ac_vals) {
  old_ls <- old_stack[[match(ac_amount, ac_vals)]]
  new_ls <- new_stack[[match(ac_amount, ac_vals)]]
  
  hist_df <- tibble(nlm_fbm = raster::getValues(old_ls),
                    fbm_fft = raster::getValues(new_ls)) %>% 
    pivot_longer(1:2, names_to = "ls_fun",
                 values_to = "value")
  
  # Overlay
  print(
    ggplot(hist_df, aes(x = value, fill = ls_fun)) +
    geom_histogram(position = "identity", alpha = .5, binwidth = 0.04) +
    labs(title = paste0("ac_amount = ", ac_amount)) +
    scale_fill_manual(values = c("midnightblue", "darkred")) +
    theme_bw()
  )
  
  # Side by side
  print(
    ggplot(hist_df, aes(x = value, fill = ls_fun)) +
      geom_histogram(position = "identity", binwidth = 0.04) +
      labs(title = paste0("ac_amount = ", ac_amount)) +
      scale_fill_manual(values = c("midnightblue", "darkred")) +
      facet_wrap(~ls_fun) +
      theme_bw()
  )
}





