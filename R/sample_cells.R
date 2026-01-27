sample_cells <- function(full_state,
                         method = c("all", "random", "chessboard"),
                         n_samples = NULL,
                         seed = NULL) {
  
  method <- match.arg(method)
  if (!method %in% c("all", "random", "chessboard")) {
    stop("method must be one of 'all', 'random' or 'chessboard'")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  grid   <- full_state$grid
  agents <- full_state$agents
  sim_id <- full_state$sim_id
  step   <- full_state$step
  grid_size <- full_state$grid_size
  fragmentation <- full_state$fragmentation
  clumped <- toroidal_clump(grid, directions = 4)
  

  # Create named vector for patches: name = patch ID, value = patch size
  patch_freq <- raster::freq(clumped, useNA = "no")
  patch_ids   <- patch_freq[, 1]
  patch_sizes <- patch_freq[, 2]
  patches <- setNames(patch_sizes, patch_ids)
  
  
  # Extract IDs of habitat cells only
  grid_vals <- raster::getValues(grid)
  habitat_cells <- which(!is.na(grid_vals))
  
  # Convert cell index to row/col
  coords <- raster::rowColFromCell(grid, habitat_cells)
  
  # Select samples
  samples <- switch(
    method,
    
    all = habitat_cells,
    
    random = {
      if (is.null(n_samples)) {
        stop("n_samples must be provided for random sampling")
      }
      sample(habitat_cells, n_samples)
    },
    
    chessboard = {
      # Every other cell based on parity of row + column
      keep <- (coords[, 1] + coords[, 2]) %% 2 == 0
      habitat_cells[keep]
    }
  )
  
  # Only present species (while all species were included in GeDo_run sample output)
  species_seq <- sort(unique(agents$species_id))
  
  # Build output
  out <- vector("list", length(samples))
  
  # Iterate through all selected samples
  for (i in seq_along(samples)) {
    
    # Cell ID
    cell <- samples[i]
    
    # Assess coordinates and patch ID
    xyloc <- raster::rowColFromCell(grid, cell)
    pid <- clumped[cell]
    
    # Abundance vector for each species
    species_counts <- sapply(
      species_seq,
      function(sp) {
        sum(
          agents$species_id == sp &
            agents$x_loc == xyloc[1] &
            agents$y_loc == xyloc[2]
        )
      }
    )
    
    # Add sample to the output list
    out[[i]] <- c(
      sim_id = sim_id,
      step = step,
      fragmentation = fragmentation,
      grid_size = grid_size,
      sample_id = i,
      loc_x = xyloc[1],
      loc_y = xyloc[2],
      patch_id = pid,
      patch_size = patches[[as.character(pid)]],
      species_counts
    )
  }
  
  # Turn list into a data frame
  out_df <- as.data.frame(do.call(rbind, out))
  colnames(out_df) <- c(
    "sim_id", "sample_id", "step", "fragmentation", "grid_size", "loc_x", "loc_y",
    paste0("sp_", species_seq)
  )
  
  return(out_df)
}
