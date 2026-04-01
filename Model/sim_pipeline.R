library(here)
library(data.table)
library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)

# Model source files
source(here("Model/parameters.R"))
source(here("Model/src/clean_run.R"))
source(here("Model/src/initialize.R"))
source(here("Model/src/landscape.R"))
source(here("Model/src/generate_agents.R"))
source(here("Model/src/distribute_agents.R"))
source(here("Model/src/record_state.R"))
source(here("Model/src/run_model_step.R"))
source(here("Model/src/birth.R"))
source(here("Model/src/death.R"))
source(here("Model/src/disperse.R"))
source(here("Model/src/immigration.R"))
source(here("Model/src/fragmentation.R"))

# Sampling and analysis source files
source(here("R/toroidal_dist.R"))
source(here("R/toroidal_clump.R"))
source(here("R/sample_cells.R"))
source(here("R/dist_decay.R"))
source(here("R/sSBR.R"))


# Multiple runs in array job ------------------------------------------------

# For each run, we need to build a var_par with varied ac and frag values,
# to get all combinations of the two parameters.

var_par_df <- tidyr::expand_grid(
  ac = seq(from = 0.1, to = 0.9, by = 0.2),
  frag = c(0.2, 0.5, 0.8),
  hab = 0.15,
  nb = 0.1,
  disp = 1,
  disp_dist = 2,
  edge = 1
)

# Parse command line argument for array index (if running on cluster)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No index provided. Please provide an index for the parameter combination to run.")
}

index <- as.numeric(args[1])
sim_id <- index
cat("Running simulation", sim_id, "with ac =", var_par_df$ac[index], "and frag =", var_par_df$frag[index], "\n")
var_par <- var_par_df[index, ]

results <- clean_run(
  mod_par = mod_par,
  var_par = var_par,
  switch = switch,
  sim_id = sim_id,
  seed = master_seed + index, # Use a different seed for each run
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
    file = paste0(here("data-raw/sampled_data/"), "sim_", sim_id, "_", 
                  "ac_", var_par$ac, 
                  "_frag_", var_par$frag, "_", 
                  "samp_", method, ".csv"),
    na = "NA"
  )
}

cat("\n")


