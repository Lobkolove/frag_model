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

#' Simulate 2D fractional Brownian motion (fBm) landscape using FFT
#'
#' Generates a fractal neutral landscape with controllable autocorrelation amount
#' using spectral synthesis and inverse FFT. The output can be returned as a matrix
#' or a spatial raster object.
#'
#' @param gr_size Integer. Size of the square grid (number of rows and columns).
#' @param resolution Numeric. Resolution of each grid cell in spatial units. Default is 1.
#' @param ac_amount Numeric. Autocorrelation amount controlling landscape smoothness.
#'   Ranges between 0 (minimum autocorrelation, highly heterogeneous) and 1 (maximum autocorrelation, smooth).
#' @param alpha.min Numeric. Minimum spectral exponent controlling heterogeneity (default 0).
#' @param alpha.max Numeric. Maximum spectral exponent controlling smoothness (default 4).
#' @param seed Optional integer. Random seed for reproducible results.
#' @param raster Logical. If TRUE, the function returns a RasterLayer object; otherwise returns a matrix.
#' @param rescale Logical. If TRUE and raster = TRUE, rescales raster values to [0, 1] range using landscapetools.
#'
#' @return Either a numeric matrix (if raster = FALSE) or a RasterLayer object (if raster = TRUE)
#'   representing the simulated landscape.
#'
#' @details
#' The function synthesizes fractal Brownian motion-like landscapes by specifying a
#' spectral density following a power law, where the spectral exponent alpha is linearly
#' interpolated from ac_amount between alpha.min and alpha.max.
#' 
#' Setting ac_amount = 0 produces very rough, noisy landscapes with little autocorrelation,
#' while ac_amount = 1 produces very smooth landscapes with strong spatial autocorrelation.
#'
#' @examples
#' # Generate a 128x128 matrix landscape with moderate autocorrelation
#' mat_landscape <- fbm_fft(gr_size = 128, ac_amount = 0.6, raster = FALSE)
#' image(mat_landscape, col = terrain.colors(100), asp = 1)
#'
#' # Generate a raster landscape with maximum autocorrelation and rescale values
#' rast_landscape <- fbm_fft(gr_size = 256, ac_amount = 1, raster = TRUE, rescale = TRUE, seed = 42)
#' plot(rast_landscape)
#'
#' @importFrom raster raster extent res
#' @importFrom landscapetools util_rescale
#' @importFrom checkmate assert_count assert_numeric assert_true assert_logical
#' @export
fbm_fft <- function(
    gr_size = 128,
    resolution = 1,
    ac_amount = 0.7,
    alpha.min = 0,
    alpha.max = 3,
    seed = NULL,
    raster = TRUE,
    rescale = TRUE
) {
  
  # Input validation
  checkmate::assert_count(gr_size, positive = TRUE)
  checkmate::assert_numeric(resolution)
  checkmate::assert_numeric(ac_amount)
  checkmate::assert_true(ac_amount >= 0)
  checkmate::assert_true(ac_amount <= 1)
  checkmate::assert_logical(rescale)  
  
  if(!is.null(seed)) set.seed(seed)
  N <- gr_size
  
  # Define alpha
  alpha <- alpha.min + ac_amount * (alpha.max - alpha.min)
  
  fx <- ifelse(0:(N - 1) <= N / 2, 0:(N - 1), 0:(N - 1) - N)
  fy <- fx
  FX <- matrix(rep(fx, each = N), nrow = N)
  FY <- matrix(rep(fy, times = N), nrow = N)
  freq <- sqrt(FX^2 + FY^2)
  freq[1, 1] <- 1 # avoid divide-by-zero at DC
  amp <- 1 / (freq^alpha)
  
  # Complex Gaussian noise
  noise <- matrix(rnorm(N * N), nrow = N) + 1i * matrix(rnorm(N * N), nrow = N)
  f_field <- amp * noise
  f_field[1, 1] <- 0 # remove mean
  
  field <- Re(fft(f_field, inverse = TRUE))
  # field <- scale(field)
  
  # Rescale to 0-1
  if(rescale) field <- scales::rescale(field)
  
  # Convert to raster and rescale
  if (raster) {
    rast <- raster::raster(field)
    raster::extent(rast) <- c(
      0,
      ncol(rast) * resolution,
      0,
      nrow(rast) * resolution
    )
    raster::res(rast) <- resolution 
    
    return(rast)
  } else {
    return(field)
  }
}


#' @noRd
#' Deprecated prototype for simulating Gaussian random fields using
#' FFT-based spectral methods.
#' Superseded by fbm_fft, which is explicitly tailored to the model,
#' trading generality for fewer options, guaranteed toroidal
#' periodicity, and more efficient usage.
grf_fft <- function(
    nx = 128, ny = 128,
    range = 20,
    sigma = 1,
    mean = 0,
    corr_fun = c("gaussian", "exponential", "matern32", "matern52", "spherical"),
    padding = NULL,        # If NULL, auto padding = 1.5 * range (rounded up)
    gradient_fun = NULL,   # A function of x and y (matrices) for adding a trend, or NULL
    gradient_strength = 1, # Scaling factor for gradient
    seed = NULL,
    periodic = FALSE       # If TRUE, disables cropping/padding and keeps fully toroidal
) {
  if(!is.null(seed)) set.seed(seed)
  corr_fun <- match.arg(corr_fun)
  
  # ------ 1. Padding and Grid Setup ------
  if(periodic) {
    nx_pad <- nx
    ny_pad <- ny
    pad_x <- pad_y <- 0
  } else {
    # Padding to suppress wraparound artifacts
    # The paper suggests at least one correlation "range" as padding
    if(is.null(padding)) padding <- ceiling(1.5 * range)
    nx_pad <- nx + 2*padding
    ny_pad <- ny + 2*padding
    pad_x <- pad_y <- padding
  }
  
  # ------ 2. Frequency Grids ------
  # For each axis, zero-based indices, centered at zero for FFT
  fx <- ifelse((0:(nx_pad-1)) <= nx_pad/2, 0:(nx_pad-1), 0:(nx_pad-1) - nx_pad)
  fy <- ifelse((0:(ny_pad-1)) <= ny_pad/2, 0:(ny_pad-1), 0:(ny_pad-1) - ny_pad)
  
  FX <- matrix(rep(fx, each=ny_pad), nrow=ny_pad)
  FY <- matrix(rep(fy, times=nx_pad), nrow=ny_pad)
  # Euclidean distance in "pixels"
  D <- sqrt(FX^2 + FY^2)
  
  # ------ 3. Correlation Spectra ------
  # Table 1 in the paper: all use relative distance d = h/range
  d <- D / range
  S <- switch(corr_fun,
              gaussian    = exp(-3*d^2),
              exponential = exp(-3*d),
              matern32    = { s <- 4.744*d; exp(-s)*(1 + s) },
              matern52    = { s <- 5.918*d; exp(-s)*(1 + s + s^2/3) },
              spherical   = { out <- 1 - d*(1.5 - 0.5*d^2); out[d>=1] <- 0; out }
  )
  
  # Fourier transform of covariance
  F_c <- fft(S)
  # Enforce non-negative spectrum for numerical stability (see eq. 19 in paper)
  F_c <- pmax(Re(F_c), 0)
  
  # ------ 4. Simulate White Noise ------
  wn <- matrix(rnorm(nx_pad*ny_pad), nrow=ny_pad)
  wn_fft <- fft(wn)
  
  # ------ 5. Apply filter in spectral space ------
  field_fft <- sqrt(F_c) * wn_fft
  field <- Re(fft(field_fft, inverse=TRUE)/(nx_pad*ny_pad))
  
  # ------ 6. Crop to remove padding ------
  if(!periodic) {
    cidx_x <- (pad_x+1):(pad_x+nx)
    cidx_y <- (pad_y+1):(pad_y+ny)
    field <- field[cidx_y, cidx_x]
  }
  
  # ------ 7. Standardize ------
  field <- mean + sigma * scale(as.vector(field))
  grf <- matrix(field, nrow=ny, ncol=nx)
  
  # ------ 8. Add gradient if desired ------
  if(!is.null(gradient_fun)) {
    xmat <- matrix(rep(seq(0,1,length.out=nx), each=ny), nrow=ny)
    ymat <- matrix(rep(seq(0,1,length.out=ny), times=nx), nrow=ny)
    grad <- gradient_fun(xmat, ymat)
    grf <- grf + gradient_strength * grad
  }
  
  return(grf)
}