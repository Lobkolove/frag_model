# Script used to generate a raster object from the old function (nlm_fbm), 
# to be compared with the output of the new function used (fbm_fft).

# Note that due to dependency issues this script can only be run on R version
# 4.4.3 or lower. Some packages used are not supported anymore and could not be
# installed from source on R version 4.5.1 or higher (compilation failed).
# Packages landscapetools v.0.5.0, RandomFieldsUtils v.1.2.5 and RandomFields
# v.3.3.14 were all downloaded from the CRAN archive and installed from source.

source("Model/src/landscape.R")

# Set parameters:
gr_size <- 128
ac_vals <- seq(0.01, 0.91, 0.10)
resolution <- 1
set.seed(42)

for (ac_amount in ac_vals) {
  # Generate grid as done in the model so far:
  env_grid <- nlm_fbm(gr_size, gr_size, 
                      resolution = resolution, 
                      fract_dim = 2 * (ac_amount))
  
  # Add information about used parameters to the filename: 
  filename <- paste0("data-raw/nlm_fbm_ac", ac_amount, ".tif")
  
  # write as GeoTIFF
  raster::writeRaster(env_grid,
                      filename = filename,
                      format   = "GTiff",
                      overwrite = TRUE)
}

