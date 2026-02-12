#' Toroidal (periodic) Euclidean distance
#'
#' Drop-in replacement for stats::dist() on a 2D (quadratic) toroidal grid.
#'
#' @param coords Numeric matrix or data.frame with two columns (x, y)
#' @param grid_size Numeric scalar giving grid side length
#'
#' @return An object of class "dist"
#' @export
toroidal_dist <- function(coords, grid_size) {
  
  coords <- as.matrix(coords)
  
  # Compute all pairwise differences
  dx <- abs(outer(coords[, 1], coords[, 1], "-"))
  dy <- abs(outer(coords[, 2], coords[, 2], "-"))
  
  # Apply toroidal wrapping by taking the minimum distance considering wrap-around
  dx <- pmin(dx, grid_size - dx)
  dy <- pmin(dy, grid_size - dy)
  
  # Compute distances
  dist_matrix <- sqrt(dx^2 + dy^2)
  
  # Extract lower triangle as dist object
  t_dist <- as.dist(dist_matrix)

  return(t_dist)

}
