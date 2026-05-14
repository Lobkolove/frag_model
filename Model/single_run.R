# Script for single local run of the model
library(here)
library(fs)
library(data.table)
library(raster)
library(dplyr)
library(tidyr)

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
source(here("R/toroidal_clump.R"))
source(here("R/sample_cells.R"))


# Simulation -------------------------------------------------------------

# Paths
output_dir <- here("output")
log_file <- here(output_dir, "simulations_log.csv")
state_dir <- here(output_dir, "model_states")
sampled_dir <- here(output_dir, "sampled_data")

# Assign a unique sim_id for this run (last used sim_id + 1)
last_sim_id <- 0
existing_log <- NULL
if (file.exists(log_file)) {
  existing_log <- fread(log_file)
  last_sim_id <- max(as.numeric(existing_log$sim_id), 0, na.rm = TRUE)
}
sim_id <- sprintf("%04d", last_sim_id + 1)

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

cat("Running simulation", sim_id, "with ac =", var_par$ac, "and frag =", var_par$frag, "\n")

# Run the model and record states at specified steps
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

# Similarly to what we did with sim_id, we can also assign a replicate number for this scenario  
scenario <- paste0("ac", var_par$ac, "_frag", var_par$frag, "_edge", var_par$edge, "_disp", var_par$disp_dist)
replicate_num <- 1
if (!is.null(existing_log)) {
  scenario_reps <- existing_log[scenario_key == scenario]
  if (nrow(scenario_reps) > 0) replicate_num <- max(scenario_reps$replicate_num) + 1
}

# We can now save all recorded steps into a single RDS file
filename <- paste0(
  sim_id,
  "_",
  scenario,
  "_r", sprintf("%03d", replicate_num)
)

state_file <- file.path(state_dir, paste0(filename, ".rds"))
saveRDS(results, state_file)
cat("✓ Simulation", sim_id, "completed and saved to:", "\n", state_file)


# Sampling ---------------------------------------------------------------

# For sampling, we want to specify which steps we want to sample from and which sampling method(s) we want to use.
# For example, here we sample from all recorded steps using random sampling only. 
recorded_steps <- names(results)
sampling_methods <- "random"
sampled_dir <- here(output_dir, "sampled_data")
sampled_files <- character()

# For each sampling method we will create a CSV file with the sampled data from all recorded steps.
for (method in sampling_methods) {

  sampled <- list()

  # Since we want to bind the sampled data from all steps into a single table, it's easier
  # to generate the single step samples in long format. 
  for (step in recorded_steps) {
    sampled[[step]] <- sample_cells(
      results[[step]],
      method = method,
      n_samples = 30, # remember to specify if using random sampling
      format = "long"
    )
  }

  sampled_long <- rbindlist(sampled)

  # Add replicate number to the sampled data.
  sampled_long[, `:=`(replicate_num = replicate_num)]

  # Here we are exporting to long format, since it is more flexible for spatial analyses.
  # Yet, other analyses may require wide format. In that case, we could reformat the data before exporting.

  # Export to CSV
  sampled_file <- file.path(sampled_dir, paste0(filename, "_samp_", method, ".csv"))
  fwrite(sampled_long, sampled_file, na = "NA")
  sampled_files <- c(sampled_files, path_rel(sampled_file, start = here()))

  cat("✓ Sampled data for method", method, "saved to:","\n", sampled_file, "\n")
}


# Logging ----------------------------------------------------------------

# Finally, we can log the simulation details into the simulations log.
# This is important for keeping track of what simulations we have run, the parameters used, 
# and where the outputs are saved. 
# It also helps to avoid overwriting existing simulations and to identify which replicates belong to which scenarios.
final_log <- data.table(
  sim_id = sim_id, 
  job_id = "local", 
  scenario_key = scenario, 
  replicate_num = replicate_num,
  run_date = Sys.Date(), 
  project_version = "1.1",
  master_seed = master_seed, 
  ac_amount = var_par$ac, 
  fragmentation = var_par$frag, 
  habitat = var_par$hab, 
  niche_breadth = var_par$nb,
  dispersal = var_par$disp, 
  dispersal_dist = var_par$disp_dist, 
  edge = var_par$edge,
  status = "complete",
  state_file = path_rel(state_file, start = here()),
  sampled_files = paste(sampled_files, collapse = "; ")
)
fwrite(final_log, log_file, append = file.exists(log_file), logical01 = TRUE)

