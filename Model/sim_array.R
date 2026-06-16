# Script to run an array job of simulations with different parameter combinations (cluster)
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

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript sim_pipeline.R <task_id> <job_id>")
}

task_id <- as.integer(args[1])
job_id <- as.numeric(args[2])

# Paths
log_file <- here("output", "simulations_log.csv")
state_dir <- here("output", "model_states")
sampled_dir <- here("output", "sampled_data")


# Simulation -------------------------------------------------------------

# Assign a unique sim_id for this run (last used sim_id + task_id)
sim_id <- unique_sim_id(log_file, increment = task_id)

# Parameters can be changed in parameters.R.
# For array jobs, it is helpful to define a grid of parameter combinations 
# and then select the combination based on the task_id.
var_par_df <- expand.grid(
  ac = 0.7,
  frag = c(0.2, 0.5, 0.8),
  hab = 0.15,
  nb = 0.1,
  disp = 1,
  disp_dist = c(1, 5, 10),
  edge = 1
)
# Make sure that the number of combinations matches the number of tasks in the array job (or a multiple thereof)!

# Determine which parameter combination to run based on task_id
n_scenarios <- nrow(var_par_df)
idx <- ((task_id - 1) %% n_scenarios) + 1
var_par <- var_par_df[idx, ]

cat(
  "Running simulation ", sim_id, " with:\n   ",
  paste(sprintf("%-10s = %s", names(var_par), unlist(var_par)),
        collapse = "\n   "),
  "\n\n", sep = ""
)

# Run the model and record states at specified steps
# master seed == last 3 digits of job_id + task_id to ensure unique seeds across runs
master_seed <- job_id %% 1000 + task_id
results <- clean_run(
  mod_par = mod_par,
  var_par = var_par,
  switch = switch,
  sim_id = sim_id,
  seed = master_seed,
  record_steps = c(
    "start",
    "pre_fragmentation",
    "post_fragmentation",
    "final"
  )
)

# Extract metadata from the recorded states to use for filenaming and logging.
# We can use the post_fragmentation step, since it is always recorded and contains all relevant metadata.
meta <- results[["post_fragmentation"]]$meta

scenario <- scenario_key(meta = meta)
replicate_num <- rep_number(meta = meta)
filename <- sim_filename(sim_id, scenario, replicate_num)

# We can now save all recorded steps into a single RDS file
state_file <- paste0(state_dir, "/", filename, ".rds")
saveRDS(results, file = state_file)
cat("\nSimulation completed.\n\nFull states were recorded for time steps [", paste(names(results), collapse = ", "), "] and saved to:\n   ", state_file, "\n", sep = "")


# Sampling ---------------------------------------------------------------

recorded_steps <- names(results)
sampled_files <- character()
sampling_methods <- list(all = "all", cb = "checkerboard", rand = "random")

for (method in names(sampling_methods)) {

  sampled <- list()

  for (step in recorded_steps) {
    sampled[[step]] <- sample_cells(
      results[[step]],
      method = sampling_methods[[method]],
      n_samples = 30,
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
# The helper function automatically writes the log entry to the specified log file.
# Since every simulation has a unique sim_id, we shouldn't have concurrency issues.
final_log <- log_entry(
  meta = meta,
  job_id = job_id, 
  scenario_key = scenario, 
  replicate_num = replicate_num,
  project_version = "1.2",
  status = "complete",
  state_file = state_file,
  sampled_files = sampled_files,
  log_file = log_file
)