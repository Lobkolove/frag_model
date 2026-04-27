chull_area <- function(
  data,
  x_col = "x_loc",
  y_col = "y_loc"
) {
  if (nrow(data) < 3) {
    return(0)
  }

  x <- data[[x_col]]
  y <- data[[y_col]]

  hull_idx <- chull(x, y)
  xh <- x[hull_idx]
  yh <- y[hull_idx]

  # Shoelace formula
  area <- 0.5 * abs(sum(xh * c(yh[-1], yh[1]) - yh * c(xh[-1], xh[1])))
  return(area)
}

toroidal_chull_area <- function(data, x_col = "x_loc", y_col = "y_loc",
                                Lx = NULL, Ly = NULL, n_repl = 3) {
  if (nrow(data) < 3) return(0)
  
  x <- data[[x_col]]
  y <- data[[y_col]]
  
  if (is.null(Lx)) Lx <- diff(range(x))
  if (is.null(Ly)) Ly <- diff(range(y))
  
  min_area <- Inf
  
  # Create replicated grid: 3x3 periods around [0,0]
  for (ix in -(n_repl %/% 2):(n_repl %/% 2)) {
    for (iy in -(n_repl %/% 2):(n_repl %/% 2)) {
      xr <- (x + ix * Lx) %% Lx  # Wrap to [0,L)
      yr <- (y + iy * Ly) %% Ly
      
      hull_idx <- chull(xr, yr)
      xh <- xr[hull_idx]
      yh <- yr[hull_idx]
      
      # Shoelace (close polygon)
      area <- 0.5 * abs(sum(xh * c(yh[-1], yh[1]) - yh * c(xh[-1], xh[1])))
      
      # Only consider if hull spans <= one period (avoids giant hulls)
      span_x <- max(xh) - min(xh)
      span_y <- max(yh) - min(yh)
      if (span_x <= Lx * 1.1 && span_y <= Ly * 1.1 && area < min_area) {
        min_area <- area
      }
    }
  }
  
  return(min_area)
}

torus_hull_area <- function(df, x_col, y_col, n) {
  pts <- df[c(x_col, y_col)]
  names(pts) <- c("x", "y")
  
  shifts <- expand.grid(dx = c(-1, 0, 1), dy = c(-1, 0, 1))
  
  poly_area <- function(x, y) {
    if (length(x) < 3) return(0)
    sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1])) / 2
  }
  
  best_area <- Inf
  
  for (k in 1:nrow(shifts)) {
    # Shift original coordinates
    x_shifted <- pts$x + shifts$dx[k] * n
    y_shifted <- pts$y + shifts$dy[k] * n
    
    # Unwrap to [0, n) window
    x <- x_shifted - n * floor(x_shifted / n)
    y <- y_shifted - n * floor(y_shifted / n)
    
    h <- chull(x, y)
    h <- c(h, h[1])
    
    area <- abs(poly_area(x[h], y[h]))
    
    if (area < best_area) {
      best_area <- area
    }
  }
  
  return(best_area)
}
