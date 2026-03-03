library(raster)
library(viridis)
library(ggplot2)
source("R/export_tiled.R")
source("Model/src/landscape.R")
source("Model/src/fragmentation.R")


# Generate and export landscape figures

  ## fbm_fft ----

# Example for Hurst = 0.5
H <- 0.5
fbm <- fbm_fft(gr_size = 400, ac_amount = H, raster = T, seed = 42)
# As raster to have consistent plotting with ls_mask output

# Plot
image(fbm, asp = 1, axes = F, xlab = "", ylab = "", 
      col = viridis(100), main = paste0("H = ", H))

# Export single image
png(paste0("pics/fbm_H", H, "2.png"), width = 2400, height = 2400)
par(mar = c(0,0,0,0))  # Remove all margins
image(fbm, asp = 1, axes = F, xlab = "", ylab = "", 
      col = viridis(100))
dev.off()

# Export tiled grid (3x3)
export_tiled(
  matrix = fbm,
  output_filename = paste0("pics/fbm_tiled_H", H, ".png")
)

# Generate sequence for H from 0 to 1 (by 0.1)
H <- seq(from = 0, to = 1, by = 0.1)
for (value in H) {
  fbm <- fbm_fft(gr_size = 400, ac_amount = value, raster = F, seed = 42)
  image(fbm, asp = 1, axes = F, col = viridis(100),
        main = paste("H =", value))
  
  # num_str <- sprintf("%+04d", as.integer(round(value * 100)))
  # filename <- paste0("pics/fbm", num_str, ".png")
  # export_tiled(fbm, output_filename = filename)
}

# GIF for H from 0 to 1 to 0 (by 0.1)
gifski::save_gif({
  for (value in H) {
    fbm <- fbm_fft(gr_size = 400, ac_amount = value, raster = T, seed = 42)
    par(mar = c(2,0,3,0))
    image(fbm, asp = 1, axes = F, xlab = "", ylab = "",
          col = viridis(100),
          main = paste("Level of Autocorrelation =", value),
          cex.main = 3)
  }
  for (value in rev(H)) {
    fbm <- fbm_fft(gr_size = 400, ac_amount = value, raster = T, seed = 42)
    par(mar = c(2,0,3,0))
    image(fbm, asp = 1, axes = F, xlab = "", ylab = "",
          col = viridis(100),
          main = paste("Level of Autocorrelation =", value),
          cex.main = 3)
  }
}, gif_file = "pics/fbm_H.gif", width = 800, height = 800, delay = 0.4)

# Fragmented landscape example
fragmentation <- 0.5
frag <- ls_mask(fbm, habitat = 0.15, fragmentation = fragmentation, seed = 43)

# Base R Plot
plot(frag, asp = 1, axes = F, 
     col = viridis(100), colNA = "grey90",
     legend = F, box = F, 
     main = paste0("Fragmentation level = ", H))
# Unfortunately gets with border when exported with png() - not sure how to fix this

# Raster to data frame for ggplot
hdf_full <- raster::as.data.frame(fbm, xy = TRUE)
hdf_frag <- raster::as.data.frame(frag, xy = TRUE)

# Plot with ggplot
gg_fbm <- ggplot(hdf_full, aes(x = x, y = y, fill = layer)) +
    geom_raster() +
    scale_fill_viridis(na.value = "grey90") +
    coord_fixed() +
    theme_void() +
    guides(fill = "none")
gg_fbm

gg_frag <- ggplot(hdf_frag, aes(x = x, y = y, fill = layer)) +
    geom_raster() +
    scale_fill_viridis(na.value = "grey90") +
    coord_fixed() +
    theme_void() +
    guides(fill = "none")
gg_frag

# Export single image
id <- 1
filename <- paste0("pics/fbm_frag_", fragmentation, "_", id, ".png")
ggsave(filename, gg_frag, 
       width = 3 * ncol(frag), height = 3 * nrow(frag), 
       units = "px", dpi = 300)


  ## grf_fft ----


# Example: matern32
grf <- grf_fft(nx = 400, ny = 400, range = 60, 
               corr_fun = "matern32", periodic = T,
               seed = 42)
image(grf, asp = 1, axes = FALSE, 
      xlab = "", ylab = "", main = "", 
      col = viridis(100))

# Export single grid
png("pics/grf_fft.png", width = 2400, height = 2400)
par(mar = c(0,0,0,0))  # Remove all margins
image(
  1:ncol(grf), 1:nrow(grf), grf,
  col = viridis(100),
  axes = FALSE, xlab = "", ylab = "", main = "",
  asp = 1
)
dev.off()

# Export tiled grid
source("R/export_tiled.R")
export_tiled(
  matrix = grf,
  output_filename = "tiled_matern32.png"
)


