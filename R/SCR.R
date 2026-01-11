library(dplyr) 

SCR <- function(model_sample) { 
  # Extract species data as P/A, dropping species with no observations 
  pa_table <- model_sample %>% 
    select(starts_with("sp")) %>% 
    select(where(~sum(.) > 0)) %>% 
    mutate(across(everything(), ~ as.integer(. > 0))) 
  
  n <- nrow(model_sample) 
  
  scr_mat <- matrix(0, n, n) 
  dist_mat <- matrix(0, n, n) 
  
  # Generate matrix with pairwise spatial distance between sites 
  coords <- model_sample %>% 
    filter(rowSums(across(starts_with("sp"))) > 0) %>% 
    select(loc_x, loc_y) pair_dist <- as.matrix(dist(coords)) 
  
  for (i in 1:n) { 
    dist_to_site <- pair_dist[i, ] 
    
    # Shuffle plots, so that tied grouping is not biased by original order. 
    new_order <- sample(1:n) 
    dist_new <- dist_to_site[new_order] 
    new_order <- new_order[order(dist_new)] 
    
    # Move focal site to the front 
    new_order <- c(i, new_order[new_order != i]) 
    comm_ordered <- pa_table[new_order, ] 
    # 1 for absence, 0 for presence 
    comm_bool <- as.data.frame((comm_ordered == 0) * 1) 
    rich <- cumprod(comm_bool) 
    scr_mat[i, ] <- as.numeric(ncol(pa_table) - rowSums(rich)) 
    dist_mat[i, ] <- dist_to_site[order(dist_to_site)] } 
  
  out <- data.frame(n = 1:nrow(scr_mat), 
                    mean_dist = colMeans(dist_mat), 
                    mean_S = colMeans(scr_mat), 
                    low_S = apply(scr_mat, 2, quantile, prob = 0.025), 
                    up_S = apply(scr_mat, 2, quantile, prob = 0.975)) 
  
  #out <- list("SCR" = scr_mat, "dist" = dist_mat) 
  return(out) 
  }

  