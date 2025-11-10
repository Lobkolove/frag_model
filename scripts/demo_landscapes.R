library(here)
library(viridis)
source(here("R", "export_tiled.R"))
source(here("R", "grf_fft.R"))
source(here("R", "fbm_fft.R"))


# Generate and export landscape figures

  ## fbm_fft ----

# Example for Hurst = 0.5
H <- 0.5
fbm <- fbm_fft(n = 400, H = H, seed = 981)

# Plot
image(fbm, asp = 1, axes = F, col = viridis(100),
      main = paste0("H = ", H))

# Export single image
png(paste0("pics/fbm_H", H, ".png"), width = 2400, height = 2400)
par(mar = c(0,0,0,0))  # Remove all margins
image(
  1:ncol(fbm), 1:nrow(fbm), fbm,
  col = viridis(100),
  axes = FALSE, xlab = "", ylab = "", main = "",
  asp = 1
)
dev.off()

# Export tiled grid (3x3)
export_tiled(
  matrix = fbm,
  output_filename = paste0("pics/fbm_tiled_H", H, ".png")
)

# Generate sequence for H from 0 to 1 (by 0.1)
H <- seq(from = 0, to = 1, by = 0.1)
for (value in H) {
  fbm <- fbm_fft(n = 400, H = value)
  image(fbm, asp = 1, axes = F, col = viridis(100),
        main = paste("H =", value))
  
  num_str <- sprintf("%+04d", as.integer(round(value * 100)))
  filename <- paste0("pics/fbm", num_str, ".png")
  export_tiled(fbm, output_filename = filename)
}


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


