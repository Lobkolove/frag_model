# code adjusted from FieldSim::midpoint (version 3.2.1)

# Generates a matrix with side length of (size 2^lsc_size + 1)
# i.e. with lsc_size = 5 the size is 2^5 + 1 = 33
# spatial autocorrelation is parameterized with the hurst factor. 
# Please note that high values mean high spatial autocorrelation
# and thus LOW fragmentation, while low values mean low spatial
# autocorrelation and thus HIGH fragmentation
# That means the hurst factor is inversely related with fragmentation per se

generate_mpd_landscape <- function(lsc_size, hurst)
{

  nblevel <- lsc_size # size of the landscape: 2^nblevel + 1
  H <- hurst          # Hurst factor in range 0 - 1
  
  Z <- matrix(0, 2, 2)
  Z[1, 2] <- rnorm(1)
  Z[2, 1] <- rnorm(1)
  Z[2, 2] <- rnorm(1) * sqrt(2^H)
  level <- 1
  
  while (level <= nblevel) {
    Y <- matrix(0, 2^(level) + 1, 2^(level) + 1)
    for (l in 1:(2^(level) + 1)) {
      for (m in 1:(2^(level) + 1)) {
        if (((m/2 - floor(m/2)) != 0) & ((l/2 - floor(l/2)) != 0)) {
          Y[m, l] <- Z[((m - 1)/2 + 1), ((l - 1)/2 + 1)]
        }
      }
    }
    varC <- (1 - 1/4 * 2^(H) - 1/8 * 2^(2 * H)) * 2^(-2 * level * H + H)
    varA <- (1 - 2^(2 * H - 2))/2^(2 * level * H)
    for (m in 1:2^(level - 1)) {
      for (l in 1:2^(level - 1)) {
        pc_x <- 2 * l
        pc_y <- 2 * m
        Y[pc_x, pc_y] <- sqrt(varC) * rnorm(1) + 
          1/4 * (Y[pc_x - 1, pc_y - 1] + Y[pc_x - 
                                             1, pc_y + 1] + Y[pc_x + 1, pc_y + 1] + 
                   Y[pc_x + 1, pc_y - 1])
        Y[pc_x + 1, pc_y] <- sqrt(varA) * rnorm(1) + 
          1/2 * (Y[pc_x + 1, pc_y - 1] + Y[pc_x + 
                                             1, pc_y + 1])
        Y[pc_x, pc_y + 1] <- sqrt(varA) * rnorm(1) + 
          1/2 * (Y[pc_x - 1, pc_y + 1] + Y[pc_x + 
                                             1, pc_y + 1])#
        if (m == 1) {
          Y[pc_x, pc_y - 1] <- sqrt(varA) * rnorm(1) + 
            1/2 * (Y[pc_x - 1, pc_y - 1] + Y[pc_x + 
                                               1, pc_y - 1])
        }
        if (l == 1) {
          Y[pc_x - 1, pc_y] <- sqrt(varA) * rnorm(1) + 
            1/2 * (Y[pc_x - 1, pc_y - 1] + Y[pc_x - 
                                               1, pc_y + 1])
        }
      }
    }
    level <- level + 1
    Z <- Y
  }
  
  # Scale to range 0 - 1
  Z01 <- (Z - min(Z))/(max(Z) - min(Z))
  
  return(Z01)

}

# Function to binarize the output of generate_mdp_landscape
# with given proportion of habitat
binarize_lsc <- function(lsc_mat, prop_hab)
{
  lsc_bin <- lsc_mat
  lsc_bin[lsc_mat < quantile(lsc_mat, 1 - prop_hab)] <- 0 
  lsc_bin[lsc_mat >= quantile(lsc_mat, 1 - prop_hab)] <- 1
  return(lsc_bin)
}


# Test the function ------

# Low hurst value -> high fragmentation per se
lsc1 <- generate_mpd_landscape(lsc_size = 5, hurst = 0.1)
require(raster)
plot(raster(lsc1))
plot(raster(binarize_lsc(lsc1, 0.1)))

# High hurst value --> low fragmentation per se
lsc2 <- generate_mpd_landscape(lsc_size = 5, hurst = 0.9)
plot(raster(lsc2))
plot(raster(binarize_lsc(lsc2, 0.5)))
