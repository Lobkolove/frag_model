dist_decay <- function(model_sample, 
                       binary = FALSE, 
                       method = "bray") { 
  # Input validation missing ! 
  
  # Extract community matrix and drop empty sites
  comm_mat <- model_sample %>%
    select(starts_with("sp")) %>%
    filter(rowSums(.) > 0)
  
  # Generate matrix with pairwise spatial distance between sites 
  coords <- model_sample %>% 
    filter(rowSums(across(sp_1:sp_1000)) > 0) %>% 
    select(loc_x, loc_y) 
  
  d <- stats::dist(coords) 
  
  # Generate Matrix with pairwise similarity between sites 
  similarity <- 1 - vegan::vegdist(comm_mat, method = method, binary = binary) 
  similarity[!is.finite(similarity)] <- NA 
  
  # df with spatial distances and respective similarities, ordered by sp. dist. 
  dat_out <- data.frame(distance = as.numeric(d), 
                        similarity = as.numeric(similarity)) 
  dat_out <- dat_out[order(dat_out$distance), ] 
  class(dat_out) <- c("dist_decay", "data.frame") 
  
  return(dat_out) 
}

plot.dist_decay <- function(x, ...) { 
  graphics::plot(similarity ~ distance, 
                 data = x, 
                 las = 1, 
                 xlab = "Distance", 
                 ylab = "Similarity", 
                 main = "Distance decay", ...) 
  dd_loess <- stats::loess(similarity ~ distance, data = x) 
  pred_sim <- stats::predict(dd_loess) 
  graphics::lines(x$distance, pred_sim, col = "red", lwd = 2) 
}


