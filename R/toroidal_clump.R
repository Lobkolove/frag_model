#' Identify habitat patches on a toroidal grid
#'
#' Extends \code{raster::clump()} to toroidal (wrap-around) landscapes. Patches
#' are first delineated as on a bounded grid; patches that touch across
#' opposing edges are then merged into one.
#'
#' @param grid A \code{RasterLayer} or a numeric matrix. Every non-\code{NA}
#'   cell is habitat, \code{NA} cells are matrix. Cell values are ignored.
#' @param directions Neighbourhood rule for connecting cells, either \code{4}
#'   (orthogonal neighbours, default) or \code{8} (including diagonals).
#'
#' @details
#' \code{raster::clump()} treats both \code{NA} and \code{0} as background.
#' Because environmental values are rescaled to [0, 1], a habitat cell can have
#' a value of exactly 0 and would be dropped, which can also split a patch it
#' bridges. To avoid this, \code{grid} is binarised internally (all non-\code{NA}
#' cells set to 1) before clumping. Consequently, a literal 0 can never be used
#' to mark background; use \code{NA} for matrix cells.
#'
#' Patches are merged across the grid borders by comparing the cells on
#' opposing edges (left-right and top-bottom) and joining the patch IDs they
#' belong to. Merged patches are assigned new consecutive IDs.
#'
#' @return The patch IDs in the same format as the input: a \code{RasterLayer}
#'   if \code{grid} was a \code{RasterLayer}, otherwise a matrix. Habitat cells
#'   hold consecutive integer patch IDs, matrix cells are \code{NA}.
#'
#' @seealso \code{\link[raster]{clump}}
#'
#' @importFrom raster raster clump as.matrix getValues values
#' @importFrom stats na.omit setNames
#'
#' @export
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

  # Binarise grid to avoid cells with habitat value 0 being treated as matrix
  rast[!is.na(rast)] <- 1

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
