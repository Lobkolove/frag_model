# Function which extends raster::clump for the use on toroidal grids.
# raster::clump assumes hard edges and no wrap-around.
# toroidal_clump merges patches which touch on opposing sides of the grid.

toroidal_clump <- function(grid, directions = 4) {

  # Convert to raster (if needed), track whether input was raster
  is_raster <- inherits(grid, "RasterLayer")

  if (is_raster) {
    rast <- grid
  } else if (is.matrix(grid)) {
    rast <- raster::raster(grid)
  } else {
    stop("Input grid must be a RasterLayer or matrix")
  }

  # Standard clumping
  clumped <- raster::clump(rast, directions = directions, gaps = FALSE)
  
  cl_mat <- raster::as.matrix(clumped)
  nrow   <- nrow(cl_mat)
  ncol   <- ncol(cl_mat)
  
  # Collect patch IDs to be merged across borders
  merge_pairs <- list()
  
  # Left–right edges
  for (r in seq_len(nrow)) {
    left  <- cl_mat[r, 1]
    right <- cl_mat[r, ncol]
    
    if (!is.na(left) && !is.na(right) && left != right) {
      merge_pairs[[length(merge_pairs) + 1]] <- c(left, right)
    }
  }
  
  # Top–bottom edges
  for (c in seq_len(ncol)) {
    top    <- cl_mat[1, c]
    bottom <- cl_mat[nrow, c]
    
    if (!is.na(top) && !is.na(bottom) && top != bottom) {
      merge_pairs[[length(merge_pairs) + 1]] <- c(top, bottom)
    }
  }
  
  # If nothing to merge, return standard clumps
  if (length(merge_pairs) == 0) {
    return(clumped)
  }
  
  merge_pairs <- unique(do.call(rbind, merge_pairs))
  
  # Build equivalence classes of patch IDs
  patch_ids <- na.omit(sort(unique(raster::getValues(clumped))))
  
  adj <- lapply(patch_ids, function(x) x)
  names(adj) <- as.character(patch_ids)
  
  for (i in seq_len(nrow(merge_pairs))) {
    a <- as.character(merge_pairs[i, 1])
    b <- as.character(merge_pairs[i, 2])
    adj[[a]] <- unique(c(adj[[a]], b))
    adj[[b]] <- unique(c(adj[[b]], a))
  }
  
  # Depth-first search to find connected components
  visited <- setNames(rep(FALSE, length(adj)), names(adj))
  comp_id <- setNames(rep(NA_integer_, length(adj)), names(adj))
  comp <- 0
  
  for (id in names(adj)) {
    if (!visited[id]) {
      comp <- comp + 1
      stack <- id
      
      while (length(stack) > 0) {
        cur <- stack[length(stack)]
        stack <- stack[-length(stack)]
        
        if (!visited[cur]) {
          visited[cur] <- TRUE
          comp_id[cur] <- comp
          # Add unvisited neighbors to stack
          for (nbr in adj[[cur]]) {
            if (!visited[nbr]) {
              stack <- c(stack, nbr)
            }
          }        
        }
      }
    }
  }
  
  # Reassign patch IDs in raster
  cl_vals <- raster::values(clumped)
  idx <- !is.na(cl_vals) & as.character(cl_vals) %in% names(comp_id)
  cl_vals[idx] <- comp_id[as.character(cl_vals[idx])]
  raster::values(clumped) <- cl_vals
  
  # Return result in original format
  if (!is_raster) {
    return(raster::as.matrix(clumped))
  } else {
    return(clumped)
  }
}
