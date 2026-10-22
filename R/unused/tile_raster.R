#' Tile a RasterLayer into an n x n mosaic
#'
#' @param r A RasterLayer (from the {raster} package).
#' @param n Integer >= 1. Number of tiles per side (default 3 -> 3x3).
#' @return A RasterLayer representing an n x n tiled version of `r`.
#'
#' @details
#' This function creates shifted copies of the input raster and mosaics them
#' together with `raster::merge()`. Shifting uses the raster's extent width/height,
#' so it works in the raster's map units.
tile_raster <- function(r, n = 3) {
  stopifnot(inherits(r, "RasterLayer"))
  stopifnot(is.numeric(n), length(n) == 1, is.finite(n), n >= 1)

  n <- as.integer(n)

  # tile size in map units (uses extent, not just ncol*res)
  ex <- raster::extent(r)
  dx <- ex@xmax - ex@xmin
  dy <- ex@ymax - ex@ymin

  # flat list of shifted rasters (critical: no nested lists)
  pieces <- lapply(0:(n - 1), function(ix) {
    lapply(0:(n - 1), function(iy) {
      raster::shift(r, dx = ix * dx, dy = iy * dy)
    })
  }) |>
    unlist(recursive = FALSE)

  do.call(raster::merge, pieces)
}
