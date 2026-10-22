# Function to create a two-dimensional fractional Brownian motion neutral landscape model.
# The function was taken from the NLMR package, a part of the rOpenSci project - https://github.com/ropensci/NLMR
# The package was developed by Marco Sciaini, Matthias Fritsch, Craig Simpkins, Cédric Scherer,and Sebastian Hanß

nlm_fbm <- function(ncol,
                    nrow,
                    resolution = 1,
                    fract_dim = 1,
                    user_seed = NULL,
                    rescale = TRUE,
                    ...) {
  # Check function arguments ----
  checkmate::assert_count(ncol, positive = TRUE)
  checkmate::assert_count(nrow, positive = TRUE)
  checkmate::assert_numeric(resolution)
  checkmate::assert_numeric(fract_dim)
  checkmate::assert_true(fract_dim > 0)
  checkmate::assert_true(fract_dim <= 2)
  checkmate::assert_logical(rescale)
  
  # specify RandomFields options ----
  RandomFields::RFoptions(cPrintlevel = 0)
  RandomFields::RFoptions(spConform = FALSE)
  RandomFields::RFoptions(...)
  
  # set RF seed ----
  RandomFields::RFoptions(seed = user_seed)
  
  # formulate and simulate fBm model
  fbm_model <- RandomFields::RMfbm(
    alpha = fract_dim)
  fbm_simu <- RandomFields::RFsimulate(fbm_model,
                                       # fBm changes x and y?
                                       y = seq.int(0, length.out = ncol),
                                       x = seq.int(0, length.out = nrow),
                                       grid = TRUE)
  
  
  # transform simulation into raster ----
  fbm_raster <- raster::raster(fbm_simu)
  
  
  # specify extent and resolution ----
  raster::extent(fbm_raster) <- c(
    0,
    ncol(fbm_raster) * resolution,
    0,
    nrow(fbm_raster) * resolution
  )
  
  # Rescale values to 0-1 ----
  if (rescale == TRUE) {
    fbm_raster <- landscapetools::util_rescale(fbm_raster)
  }
  
  return(fbm_raster)
}