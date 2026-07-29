ls_mask <- function(grid,
                    habitat,
                    fragmentation,
                    seed = NULL) {
  
  # Input validation
  is_raster <- inherits(grid, "RasterLayer")
  if (!is_raster) {
    if (is.matrix(grid)) {
      grid <- raster::raster(grid)
    } else {
      stop("Input grid must be a RasterLayer or matrix")
    }
  }
  if (grid@ncols != grid@nrows) stop("Grids with different numbers of rows and columns are not supported yet")
  if (habitat <= 0 || habitat >= 1) stop("'habitat' must be a value between 0 and 1")
  if (fragmentation < 0 || fragmentation > 1) stop("'fragmentation' must be a value between 0 and 1")
  
  # Generate full grid which will be turned into mask
  mask_prep <- fbm_fft(gr_size = grid@ncols, 
                       ac_amount = 1 - fragmentation, 
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
  
  # Apply mask to input grid
  fragmented_grid <- raster::mask(grid, mask, maskvalue = 1)
  
  if (is_raster) {
    return(fragmented_grid)
  } else {
    return(raster::as.matrix(fragmented_grid))
  }
}


fragment <- function(full_state,
                     habitat,
                     fragmentation,
                     seed = NULL,
                     ...) {
  
  grid <- full_state$grid
  agents <- full_state$agents
  ss_abund <- full_state$ss_abund
  
  # Apply fragmentation to grid
  fragmented_grid <- ls_mask(grid = grid,
                             habitat = habitat,
                             fragmentation = fragmentation,
                             seed = seed)
  
  # Assess which cells are habitat and edge cells (useful for faster processes later)
  if (inherits(fragmented_grid, "RasterLayer")) {
    habitat_cells <- raster::Which(!is.na(fragmented_grid), cells = TRUE)
    edges <- boundaries(fragmented_grid, type = "inner", asNA = TRUE)
  } else {
    habitat_cells <- which(!is.na(fragmented_grid))
    edges <- boundaries(raster::raster(fragmented_grid), type = "inner", asNA = TRUE)
  }
  
  # Only keep agents which are on habitat cells 
  keep <- !is.na(fragmented_grid[cbind(agents$x_loc, agents$y_loc)])
  agents <- agents[keep, , drop = FALSE]

  # Update ss_abund (turn matrix cells to NA)
  ss_abund <- ss_abund %>%
    dplyr::mutate(
      habitat = !is.na(fragmented_grid[cbind(x_loc, y_loc)]),
      n = ifelse(habitat, n, NA)
    ) %>%
    dplyr::select(-habitat)

  meta <- full_state$meta
  meta$fragmentation <- fragmentation
  meta$habitat <- habitat
  
  record_state(
    model_state = record_step(
      grid = fragmented_grid,
      agents = agents,
      step = full_state$step,
      cells = list(
        habitat = habitat_cells,
        edges = edges
      )
    ),
    meta = meta,
    step_label = "post_fragmentation"
  )
}

