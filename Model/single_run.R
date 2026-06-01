# Script for single local run of the model
library(here)
library(fs)
library(data.table)
library(raster)
library(dplyr)
library(tidyr)
library(igraph)

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
source(here("Model/src/registry.R"))

# Sampling source files
source(here("R/toroidal_clump.R"))
source(here("R/sample_cells.R"))

# Paths
log_file <- "output/simulations_log.csv"
state_dir <- "output/model_states"
sampled_dir <- "output/sampled_data"


# Simulation -------------------------------------------------------------

# Assign a unique sim_id for this run (last used sim_id + 1)
sim_id <- unique_sim_id(log_file)

# Parameters can be changed in parameters.R.
# If needed, we can change them here for testing purposes. For example:
master_seed <- 12345
var_par <- list(
  ac = 0.7,
  frag = 0.5,
  hab = 0.15,
  nb = 0.1,
  disp = 1,
  disp_dist = 2,
  edge = 1
)

cat(
  "Running simulation ", sim_id, " with:\n   ",
  paste(sprintf("%-10s = %s", names(var_par), unlist(var_par)),
        collapse = "\n   "),
  "\n\n", sep = ""
)

# Run the model and record states at specified steps
# No need to specify master_seed here, since it is already defined in parameters.R
results <- clean_run(
  mod_par = mod_par,
  var_par = var_par,
  switch = switch,
  sim_id = sim_id,
  seed = master_seed,
  record_steps = c(
    # "start",
    # "pre_fragmentation",
    "post_fragmentation",
    "final"
  )
)

# Extract metadata from the recorded states to use for filenaming and logging.
# We can use the post_fragmentation step, since it is always recorded and contains all relevant metadata.
meta <- results[["post_fragmentation"]]$meta

# Similarly to what we did with sim_id, we can also assign a replicate number for this scenario  
scenario <- scenario_key(meta = meta)
replicate_num <- rep_number(meta = meta)
filename <- sim_filename(sim_id, scenario, replicate_num)

# We can now save all recorded steps into a single RDS file
state_file <- paste0(state_dir, "/", filename, ".rds")
saveRDS(results, state_file)
cat("\nSimulation completed.\n\nFull states were recorded for time steps [", paste(names(results), collapse = ", "), "] and saved to:\n   ", state_file, "\n", sep = "")


# Sampling ---------------------------------------------------------------

# For sampling, we want to specify which steps we want to sample from and which sampling method(s) we want to use.
# For example, here we sample from all recorded steps using random sampling only. 
recorded_steps <- names(results)
sampled_files <- character()
sampling_methods <- list(rand = "random")

# For each sampling method we will create a CSV file with the sampled data from all recorded steps.
for (method in names(sampling_methods)) {

  sampled <- list()

  # Since we want to bind the sampled data from all steps into a single table, it's easier
  # to generate the single step samples in long format. 
  for (step in recorded_steps) {
    sampled[[step]] <- sample_cells(
      results[[step]],
      method = sampling_methods[[method]],
      n_samples = 30, # remember to specify if using random sampling
      format = "long"
    )
  }

  sampled_long <- rbindlist(sampled)

  # Here we are exporting to long format, since it is more flexible for spatial analyses.
  # Yet, other analyses may require wide format. In that case, we could reformat the data before exporting.

  # Export to CSV
  sampled_file <- paste0(sampled_dir, "/", filename, "_", method, ".csv")
  fwrite(sampled_long, sampled_file, na = "NA")
  sampled_files <- c(sampled_files, sampled_file)

}

cat("\nData was sampled for all time steps using methods [", paste(sampling_methods, collapse = ", "), "] and saved to:\n   ", paste(sampled_files, collapse = "\n   "), "\n", sep = "")

# Logging ----------------------------------------------------------------

# Finally, we can log the simulation details into the simulations log.
# This is important for keeping track of what simulations we have run, the parameters used, 
# and where the outputs are saved. 
# It also helps to avoid overwriting existing simulations and to identify which replicates belong to which scenarios.
log_entry(
  meta = meta,
  job_id = "local", 
  scenario_key = scenario, 
  replicate_num = replicate_num,
  project_version = "1.1",
  status = "complete",
  state_file = state_file,
  sampled_files = sampled_files
)

