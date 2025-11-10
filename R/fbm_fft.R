#' @noRd
fbm_fft <- function(
    n = 128, 
    H = 0.7, 
    alpha.min = 0, alpha.max = 4, 
    seed = NULL
    ) {
  if(!is.null(seed)) set.seed(seed)
  alpha <- alpha.min + H * (alpha.max - alpha.min)
  nx <- n; ny <- n
  
  # Frequency grid
  fx <- ifelse(0:(nx-1) <= nx/2, 0:(nx-1), 0:(nx-1)-nx)
  fy <- ifelse(0:(ny-1) <= ny/2, 0:(ny-1), 0:(ny-1)-ny)
  FX <- matrix(rep(fx, each=ny), ny)
  FY <- matrix(rep(fy, times=nx), ny)
  freq <- sqrt(FX^2 + FY^2)
  freq[1,1] <- 1 # prevent div by zero at DC
  amp <- 1/(freq^alpha)
  
  # Gaussian noise in complex FFT space
  noise <- matrix(rnorm(nx*ny), nrow=ny) + 1i * matrix(rnorm(nx*ny), nrow=ny)
  f_field <- amp * noise
  f_field[1,1] <- 0 # Optionally: zero mean
  
  field <- Re(fft(f_field, inverse=TRUE))
  field <- scale(field) # Standardize to mean 0, sd 1
  return(field)
}

