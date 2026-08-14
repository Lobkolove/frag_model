expand_hull_area <- function(
  data, 
  S, 
  target_area,
  x_col = "x_loc",
  y_col = "y_col",
  n_iter = 500
) {
  N <- nrow(data)
  A_cur <- chull_area(data[S, ], x_col, y_col)
  E_cur <- abs(A_cur - target_area)
  
  best_S <- S
  best_E <- E_cur
  best_A <- A_cur
  
  for (i in seq_len(n_iter)) {
    out_idx <- sample(S, 1)
    in_cand <- setdiff(seq_len(N), S)
    if (length(in_cand) == 0) break
    in_idx <- sample(in_cand, 1)
    
    S_new <- setdiff(S, out_idx)
    S_new <- c(S_new, in_idx)
    
    A_new <- chull_area(data[S_new, ], x_col, y_col)
    E_new <- abs(A_new - target_area)
    
    # greedy: accept if closer to target
    if (E_new < E_cur) {
      S <- S_new
      A_cur <- A_new
      E_cur <- E_new
      
      if (E_new < best_E) {
        best_E <- E_new
        best_S <- S_new
        best_A <- A_new
      }
    }
  }
  
  best_S
}

chull_area_sample <- function(
  data, 
  target_area,
  n_sample = 30,
  x_col = "x_loc",
  y_col = "y_loc",
  A_max = target_area * 1.2,
  A_min = target_area * 0.8,
  max_tries = 100,
  n_iter = 500
) {

  N <- nrow(data)
  if (N < n_sample) stop("Not enough habitat cells to sample.")
    
  for (try_i in seq_len(max_tries)) {
    # start from a random focal cell
    S <- sample(N, 1)
    
    while (length(S) < n_sample) {
      # candidates: all cells not yet in S
      cand <- setdiff(seq_len(N), S)
      if (length(cand) == 0) break
      
      # try adding a random candidate
      new_idx <- sample(cand, 1)
      S_try <- c(S, new_idx)
      
      A_try <- chull_area(data[S_try, ], x_col, y_col)
      
      # accept only if area stays below A_max
      if (A_try <= A_max) {
        S <- S_try
      }
      # else: reject this candidate and try another
    }
    
    if (length(S) == n_sample) {
      
      A_final <- chull_area(data[S, ], x_col, y_col)
      
      # If too small, try to expand toward target
      if (A_final < target_area) {
        S <- expand_hull_area(data, S, target_area, x_col, y_col, n_iter)
        A_final <- chull_area(data[S, ], x_col, y_col)
      }
      
      # Accept if reasonably close
      if (A_final >= A_min && A_final <= A_max) {
        return(data[S, ])
      }
    }
  }
  
  warning("Could not find a sample with area within the specified range after ", max_tries, " tries. Returning the best found sample.")
  
  return(data[S, ])
}

