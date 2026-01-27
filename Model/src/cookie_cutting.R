cookie_cutting <- function(grid = grid,
                           agents = agents,
                           agents_grid = agents_grid,
                           habitat = mod_par$habitat_percent,
                           fragmentation = mod_par$frag_factor) {

  # frag_grid <- nlm_fbm(mod_par$grid_size, mod_par$grid_size, resolution = 1, fract_dim = 2 * (1 - fragmentation))
  frag_grid <- fbm_fft(gr_size = mod_par$grid_size, 
                       ac_amount = (1 - fragmentation), 
                       raster = T, seed = seed)

  suppressWarnings(binary_grid <- landscapetools::util_binarize(frag_grid, habitat)) # Binarize fragmentation grid as only matrix/habitat values are needed
  binary_grid[binary_grid == 1] <- NA # Subset matrix to NA


  binary_grid[binary_grid > 1] <- grid[binary_grid > 1] # Merge the layers into final space grid
  binary_grid[binary_grid == 0] <- 0.001

  delete_agents <- vector()

  clumped <- raster::clump(binary_grid, directions = 4)
  patch_matrix <- as.matrix(clumped)
  for (i in 1:nrow(agents)) {
    cur_loc <- c(agents$x_loc[i], agents$y_loc[i])
    if (is.na(binary_grid[cur_loc[1], cur_loc[2]])) {
      delete_agents <- append(delete_agents, i)
    }
  }
  agents <- agents[-c(delete_agents), ]

  # assign patch number to agent

  for (j in 1:nrow(agents)) {
    agents$patch_id[j] <- patch_matrix[agents$x_loc[j], agents$y_loc[j]]
  }
  agents_grid[is.na(binary_grid)] <- NA

  return_list <- list(grid = binary_grid, agents = agents, agents_grid = agents_grid)
  return(return_list)
}