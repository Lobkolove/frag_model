sim_id <- 1
ac_levels <- seq(0.1, 0.9, by = 0.2)
frag_levels <- c(0.2, 0.5, 0.8)

data_list <- list()

# For each ac level, read in all sampled data and merge
for (ac in ac_levels) {
  # Initialize an empty list to store data frames
  tmp_list <- list()
  
  # Loop through each frag level and read the corresponding CSV file
  for (frag in frag_levels) {
    file_name <- paste0("data-raw/sampled_data/sim_", sim_id, "_ac_", ac, "_frag_", frag, "_samp_random.csv")
    tmp_list[[length(tmp_list) + 1]] <- read.csv(file_name)
  }
  
  # Combine all data frames into one
  combined_data <- do.call(rbind, tmp_list)

  # Add the combined data frame to the main list
  data_list[[as.character(ac)]] <- combined_data
}


