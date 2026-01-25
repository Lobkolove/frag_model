######################## Grid Generating Function ##############################

# The function takes 4 parameters from the 'par' data.frame and generates 2
# landscapes using a fractional Brownian motion model.
# 1 grid simulates an auto correlated landscape and the other simulates
# fragmented landscape.
# The 2 grids are merged to create the final simulation space. The function
# returns a raster object.

generate_grid <- function(gr_size, ac_amount, frag_amount, hab_amount) {
  # if (ac_amount < 0.8) { # Generate auto correlated landscape (if else statement is for cases where fractal dimension value is high according to nlm_fbm documentation)
  #   env_grid <- nlm_fbm(gr_size, gr_size, resolution = 1, fract_dim = 2 * (ac_amount)) # multiply by 2 to simplify parameter settings (from 0-1 to 0-2)
  # } else {
  #   env_grid <- nlm_fbm(gr_size, gr_size, resolution = 1, fract_dim = 2 * (ac_amount), modus_operandi = "easygoing")
  # }
  # if (frag_amount > 0.2) { # Generate fragmented landscape (if else statement is for cases where fractal dimension value is high according to nlm_fbm documentation)
  #   frag_grid <- nlm_fbm(gr_size, gr_size, resolution = 1, fract_dim = 2 * (1 - frag_amount)) # Multiply and invert value to simplify parameter settings
  # } else {
  #   frag_grid <- nlm_fbm(gr_size, gr_size, resolution = 1, fract_dim = 2 * (1 - frag_amount), modus_operandi = "easygoing")
  # }
    # env_grid <- nlm_fbm(gr_size, gr_size, resolution = 1, fract_dim = 2 * (ac_amount))
  # frag_grid <- nlm_fbm(gr_size, gr_size, resolution = 1, fract_dim = 2 * (1 - frag_amount))
  # 
  # if (hab_amount < 1) {
  #   suppressWarnings(binary_grid <- landscapetools::util_binarize(frag_grid, hab_amount)) # Binarize fragmentation grid as only matrix/habitat values are needed
  #   binary_grid[binary_grid == 1] <- NA # Subset matrix to NA
  # 
  # 
  #   binary_grid[binary_grid > 1] <- env_grid[binary_grid > 1] # Merge the layers into final space grid
  #   binary_grid[binary_grid == 0] <- 0.001
  # } else {
    # binary_grid <- env_grid
  # }
  # return(binary_grid)
  return(fbm_fft(gr_size = gr_size, ac_amount = ac_amount, raster = T))
}