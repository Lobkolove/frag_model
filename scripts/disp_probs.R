source("Model/src/disperse.R")

grid_size <- 50

n_sims <- 1e6

sd_dist <- 2
mean_dist <- c(1, 2, 4, 8)

for (m in mean_dist) {

  count <- 0
  for (i in 1:n_sims) {
    cur_loc <- c(sample(1:grid_size, 1), sample(1:grid_size, 1))
    new_loc <- toroidal_disperse(
      cur_loc, 
      d_sd = sd_dist, 
      d_mean = m,
      grid_size = grid_size,
      kernel = "exponential"
    )

    if (new_loc[1] == cur_loc[1] && new_loc[2] == cur_loc[2]) {
      count <- count + 1
    }
  }
  cat("Mean distance:", m, 
  "\nProbability of no movement:", count / n_sims * 100, "%\n\n")
}
  
# Mean distance: 1 
# Probability of no movement: 42.8231 %

# Mean distance: 2 
# Probability of no movement: 24.4157 %

# Mean distance: 4 
# Probability of no movement: 13.1358 %

# Mean distance: 8 
# Probability of no movement: 6.7788 %