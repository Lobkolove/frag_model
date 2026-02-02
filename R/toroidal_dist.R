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
  n <- nrow(coords)
  x <- coords[, 1]
  y <- coords[, 2]
  
  dx <- abs(outer(x, x, "-"))
  dx <- pmin(dx, grid_size - dx)
  
  dy <- abs(outer(y, y, "-"))
  dy <- pmin(dy, grid_size - dy)
  
  dmat <- sqrt(dx^2 + dy^2)
  
  # extract upper triangle in dist() order
  d <- dmat[upper.tri(dmat)]
  
  structure(
    d,
    Size = n,
    Labels = rownames(coords),
    Diag = FALSE,
    Upper = FALSE,
    method = "toroidal",
    call = match.call(),
    class = "dist"
  )
}
