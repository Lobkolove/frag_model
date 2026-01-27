ls_mask <- function(grid,
                    habitat,
                    fragmentation,
                    seed = NULL) {
  
  # Input validation
  if (!inherits(grid, "RasterLayer")) stop("Input grid must be a RasterLayer")
  if (grid@ncols != grid@nrows) stop("Grids with different numbers of rows and columns are not supported yet")
  if (habitat <= 0 || habitat >= 1) stop("'habitat' must be a value between 0 and 1")
  if (fragmentation < 0 || fragmentation > 1) stop("'fragmentation' must be a value between 0 and 1")
  
  # Generate full grid which will be turned into mask
  mask_prep <- fbm_fft(gr_size = grid@ncols, 
                       ac_amount = fragmentation, 
                       raster = T, 
                       seed = seed)
  
  # Force mask_prep to share extent/res/CRS with input grid
  raster::extent(mask_prep) <- raster::extent(grid)
  raster::res(mask_prep)    <- raster::res(grid)
  raster::crs(mask_prep)    <- raster::crs(grid)
  
  # Assess threshold for binarization, based on aimed habitat proportion
  threshold <- raster::quantile(mask_prep, probs = 1 - habitat)[[1]]
  
  # Turn full grid into binary mask, using this threshold
  mask <- raster::cut(mask_prep, breaks = c(-Inf, threshold, Inf))
  
  # Old version with landscapetools dependency, replaced by last 2 steps:
  # suppressWarnings(mask <- landscapetools::util_binarize(mask_prep, habitat)) 
  
  # Apply mask to input grid
  fragmented_grid <- raster::mask(grid, mask, maskvalue = 1)
  
  fragmented_grid
}


fragment <- function(grid,
                     agents,
                     agents_grid,
                     habitat,
                     fragmentation,
                     seed = NULL,
                     ...) {
  
  # Apply fragmentation to grid
  fragmented_grid <- ls_mask(grid = grid,
                             habitat = habitat,
                             fragmentation = fragmentation,
                             seed = seed)
  
  # Only keep agents which are on habitat cells 
  keep <- !is.na(fragmented_grid[cbind(agents$x_loc, agents$y_loc)])
  agents <- agents[keep, , drop = FALSE]
  
  # Assign patch number (id) to each agent
  clumped <- raster::clump(fragmented_grid, directions = 4)
  patch_matrix <- as.matrix(clumped)
  agents$patch_id <- patch_matrix[cbind(agents$x_loc, agents$y_loc)]
  
  # Update agents_grid (turn matrix cells to NA)
  agents_grid[is.na(fragmented_grid)] <- NA
  
  return_list <- list(grid = fragmented_grid,
                      agents = agents,
                      agents_grid = agents_grid)
  return(return_list)
}
