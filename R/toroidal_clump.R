# Function to use instead of raster::clump if working on a toroidal grid.
# raster::clump assumes hard edges and no wrap-around.
# toroidal_clump merges patches which touch on opposing sides of the grid.

toroidal_clump <- function(grid, directions = 4, reindex = TRUE) {
  
  # Standard clumping
  clumped <- raster::clump(grid, directions = directions, gaps = FALSE)
  
  cl_mat <- as.matrix(clumped)
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
  patch_ids <- sort(unique(merge_pairs))
  patch_ids <- patch_ids[!is.na(patch_ids)]
  
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
        cur <- stack[1]
        stack <- stack[-1]
        
        if (!visited[cur]) {
          visited[cur] <- TRUE
          comp_id[cur] <- comp
          stack <- c(stack, adj[[cur]])
        }
      }
    }
  }
  
  # Reindex (to be added yet!) and reassign patch IDs in raster
  cl_vals <- raster::values(clumped)
  idx <- !is.na(cl_vals) & as.character(cl_vals) %in% names(comp_id)
  cl_vals[idx] <- comp_id[as.character(cl_vals[idx])]
  raster::values(clumped) <- cl_vals
  
  clumped
}
