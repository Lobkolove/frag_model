library(data.table)
library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)

# Model source files
source("Model/parameters.R")
source("Model/src/clean_run.R")
source("Model/src/initialize.R")
source("Model/src/landscape.R")
source("Model/src/generate_agents.R")
source("Model/src/distribute_agents.R")
source("Model/src/record_state.R")
source("Model/src/run_model_step.R")
source("Model/src/birth.R")
source("Model/src/death.R")
source("Model/src/disperse.R")
source("Model/src/immigration.R")
source("Model/src/fragmentation.R")

# Sampling and analysis source files
source("R/toroidal_dist.R")
source("R/toroidal_clump.R")
source("R/sample_cells.R")
source("R/dist_decay.R")
source("R/sSBR.R")

sim_id <- "test_run"

# Single run with static parameters for testing
results <- clean_run(
  mod_par = mod_par,
  var_par = var_par,
  switch = switch,
  sim_id = sim_id,
  seed = master_seed,
  record_steps = c("start", "pre_fragmentation", "post_fragmentation", "final")
)

# Sample test data and merge into a single data frame for each sampling method
recorded_steps <- names(results)
sampling_methods <- c("random", "chessboard", "all")
sampled_data <- list()
for (method in sampling_methods) {
  for (step in recorded_steps) {
    sampled <- sample_cells(
      results[[step]],
      method = method,
      n_samples = 30,
      format = "long"
    )
    sampled_data[[method]][[step]] <- sampled
    print(head(sampled))
  }
}

# Combine sampled data for each method into a single df with wide format
for (method in sampling_methods) {

  sampled_data[[method]] <- rbindlist(sampled_data[[method]])
  
  sampled_data[[method]] <- sampled_data[[method]] %>%
    tidyr::pivot_wider(names_from  = species_id,
                       values_from = n,
                       values_fill = 0,
                       names_prefix = "sp_",
                       names_sort = TRUE)
  
  # Some samples had no observations of any species, resulting in NA values for species_id and count.
  # When pivoting to wide format, this can result in a column "sp_NA" which can be removed.
  if ("sp_NA" %in% colnames(sampled_data[[method]])) {
    sampled_data[[method]] <- sampled_data[[method]] %>%
      select(-sp_NA)
  }
}

# Export random samling data to CSV
data.table::fwrite(
  sampled_data[["random"]], 
  file = paste0("../data-raw/sampled_data/sim_", sim_id, "_random.csv"),
  na = "NA"
)

# Multiple runs in a loop ------------------------------------------------

# For each run, we need to build a var_par with varied ac and frag values,
# to get all combinations of the two parameters.

var_par_df <- tidyr::expand_grid(
  ac = seq(from = 0.1, to = 0.9, by = 0.1),
  frag = c(0.2, 0.5, 0.8),
  hab = 0.15,
  nb = 0.1,
  disp = 1,
  disp_dist = 2,
  edge = 1
)

for (i in seq_len(nrow(var_par_df))) {
  
  sim_id <- i
  cat("Running simulation ", sim_id, " with ac = ", var_par_df$ac[i], " and frag = ", var_par_df$frag[i], "\n")
  var_par <- var_par_df[i, ]

  results <- clean_run(
    mod_par = mod_par,
    var_par = var_par,
    switch = switch,
    sim_id = sim_id,
    seed = master_seed + i, # Use a different seed for each run
    record_steps = c("start", "pre_fragmentation", "post_fragmentation", "final")
  )

  recorded_steps <- names(results)
  sampling_methods <- "random"
  for (method in sampling_methods) {
    sampled <- list()
    for (step in recorded_steps) {
      sampled[[step]] <- sample_cells(
        results[[step]],
        method = method,
        n_samples = 30,
        format = "long"
      )
    }
    
    sampled_df <- rbindlist(sampled)
    
    sampled_df <- sampled_df %>%
      tidyr::pivot_wider(names_from  = species_id,
                         values_from = n,
                         values_fill = 0,
                         names_prefix = "sp_",
                         names_sort = TRUE)
    
    if ("sp_NA" %in% colnames(sampled_df)) {
      sampled_df <- sampled_df %>%
        select(-sp_NA)
    }
    
    # Export to CSV
    data.table::fwrite(
      sampled_df, 
      file = paste0("../data-raw/sampled_data/sim_", sim_id, "_", 
                    "ac_", var_par$ac, 
                    "_frag_", var_par$frag, "_", 
                    "samp_", method, ".csv"),
      na = "NA"
    )
  }

  cat("\n")
}
