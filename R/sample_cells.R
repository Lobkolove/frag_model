sample_cells <- function(full_state,
                         method = c("all", "random", "chessboard"),
                         n_samples = NULL,
                         format = c("wide", "long"),
                         seed = NULL) {
  
  method <- rlang::arg_match(method)
  format <- rlang::arg_match(format)
  
  if (!is.null(seed)) set.seed(seed)
  
  # Extract grid and species abundances per cell from state
  grid            <- full_state$grid
  ss_abund        <- full_state$ss_abund
  
  # Extract IDs and coordinates of habitat cells only
  grid_vals <- raster::getValues(grid)
  habitat_cells <- which(!is.na(grid_vals))
  # Convert cell index to row/col
  coords <- raster::rowColFromCell(grid, habitat_cells)
  
  # Select cells to be sampled
  sampled_cells <- switch(
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
  
  # Assess coordinates of sample cells
  sample_coords <- raster::rowColFromCell(grid, sampled_cells)
  
  samples <- data.frame(
    sample_id = seq_len(nrow(sample_coords)),
    cell_id = sampled_cells,
    x_loc = sample_coords[, 1],
    y_loc = sample_coords[, 2]
  )
  
  # Assess patch information
  clumped <- toroidal_clump(grid, directions = 4) 
  patch_freq <- raster::freq(clumped, useNA = "no")
  patch_ids <- patch_freq[, 1]
  patch_sizes <- patch_freq[, 2]
  # Create named vector for patches: name = patch ID, value = patch size
  patches <- setNames(patch_sizes, patch_ids)
  # Add patch info to samples df
  samples$patch_id <- clumped[samples$cell_id]
  samples$patch_size <- unname(patches[as.character(samples$patch_id)])
  
  # Add static metadata columns 
  samples <- samples %>% 
    dplyr::mutate(sim_id         = full_state$sim_id,
                  master_seed    = full_state$master_seed,
                  step           = full_state$step,
                  grid_size      = full_state$grid_size,
                  fragmentation  = full_state$fragmentation)
  
  
  # Merge samples df and species abundances per cell
  out_long <- dplyr::left_join(samples, 
                               ss_abund, 
                               by = c("x_loc", "y_loc"))
  
  # Return long format
  if (format == "long") return(out_long)

  # If needed, reformat and return wide format
  out_wide <- out_long %>% 
    tidyr::pivot_wider(names_from  = species_id,
                       values_from = n,
                       values_fill = 0,
                       names_prefix = "sp_",
                       names_sort = TRUE)
  return(out_wide)

}
