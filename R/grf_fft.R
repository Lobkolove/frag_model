#' @noRd
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
