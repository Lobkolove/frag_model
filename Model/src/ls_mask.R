ls_mask <- function(grid = grid,
                    agents = agents,
                    agents_grid = agents_grid,
                    habitat = mod_par$habitat_percent,
                    fragmentation = mod_par$frag_factor) {
  
  # Generate full grid which will be turned into mask
  mask_prep <- fbm_fft(mod_par$grid_size, 
                       ac_amount = fragmentation, 
                       raster = T, seed = seed)
  
  # Assess threshold for binarization, based on aimed habitat proportion
  threshold <- raster::quantile(mask_prep, probs = 1 - habitat)[[1]]
  
  # Turn full grid into binary mask, using this threshold
  mask <- raster::cut(mask_prep, breaks = c(-Inf, threshold, Inf))
  
  # Old version with landscapetools dependency, replaced by last 2 steps:
  # suppressWarnings(mask <- landscapetools::util_binarize(mask_prep, habitat)) 
  
  # Apply mask to input grid
  fragmented_grid <- raster::mask(grid, mask, maskvalue = 1)
  
  delete_agents <- vector()
  
  clumped <- raster::clump(fragmented_grid, directions = 4)
  patch_matrix <- as.matrix(clumped)
  for (i in 1:nrow(agents)) {
    cur_loc <- c(agents$x_loc[i], agents$y_loc[i])
    if (is.na(fragmented_grid[cur_loc[1], cur_loc[2]])) {
      delete_agents <- append(delete_agents, i)
    }
  }
  agents <- agents[-c(delete_agents), ]
  
  # assign patch number to agent
  
  for (j in 1:nrow(agents)) {
    agents$patch_id[j] <- patch_matrix[agents$x_loc[j], agents$y_loc[j]]
  }
  agents_grid[is.na(fragmented_grid)] <- NA
  
  return_list <- list(grid = fragmented_grid,
                      agents = agents,
                      agents_grid = agents_grid)
  return(return_list)
}