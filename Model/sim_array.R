# Script for a simulation array (usually cluster)
library(here)
library(fs)
library(data.table)
library(raster)
library(dplyr)
library(tidyr)
library(igraph)
library(collapse)
library(checkmate)
library(withr)
library(scales)

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

# Settings ---------------------------------------------------------------

# Select the time steps to be recorded and sampled.
steps_to_record <- c(
  "start",
  "pre_fragmentation",
  "post_fragmentation",
  "final"
)
# Note that if steps_to_record is empty, simulations will still run but
# no states will be recorded and no sampling can be performed.

# Should full states be exported as RDS files? 
# If yes, set export_states = TRUE.
export_states <- TRUE

# Should the results be sampled? 
# If yes, specify the sampling method(s) by uncommenting the desired method(s) below.
# If no, assign sampling_methods as empty list().
sampling_methods <- list(
  all = "all", 
  cb = "checkerboard", 
  rand = "random"
)

# Should the simulations be logged? 
# If yes, set to_log = TRUE. 
to_log <- TRUE
# All simulations relevant for analyses should be logged!
# Set to_log = FALSE only for testing and debugging purposes!

# Notes (any additional information to be printed to the console output)
# (will be included in the logs/*.out files)
notes <- NULL

# Is the current array running in parallel to another? 
# In that case we need to increment the sim_id by the number of simulations in the other array to avoid overlapping sim_ids.
# If not, set spacer = 0. If yes, set spacer = number of simulations in the other array.
spacer <- 27


# Simulation -------------------------------------------------------------

# Define a unique identifier for this simulation run (sim_id)
if (to_log) {
  # Last used sim_id is determined from the log file, and the next sim_id is assigned.
  sim_id <- unique_sim_id(log_file = log_file, increment = task_id + spacer)
} else {
  # If not logging, we can just use the task_id as sim_id
  sim_id <- task_id
}

# To avoid overlapping seeds, for arrays it's best to override the master_seed from parameters.R
# and use the job_id as root instead and increment it by the task_id. 
# This way, the chance of overlapping seeds is very low, 
# even when the number of arrays and simulations gets high.
master_seed <- job_id %% 1000 + task_id
# master seed == last 3 digits of job_id + task_id 

# Other parameters should be changed in parameters.R.
# Make sure that the number of tasks in the array job matches the number of
# parameter combinations in var_par (or a multiple thereof)! 

# Determine which parameter combination to run based on task_id
n_scenarios <- nrow(var_par)
idx <- ((task_id - 1) %% n_scenarios) + 1
cur_var_par <- var_par[idx, ]

cat(
  "Running simulation ", sim_id, " with:\n\n   ",
  "Dispersal mode: ", ifelse(switch$dispersal_type == 1, "short_long", "random"), "\n\n   ",
  "Parameters:\n      ",
  paste(sprintf("%-10s = %s", names(cur_var_par), unlist(cur_var_par)), collapse = "\n      "),
  "\n\n", 
  sep = ""
)

if (!is.null(notes)) cat("   Notes:   ", notes, "\n\n", sep = "")

# Run the model and record states at specified steps
results <- clean_run(
  mod_par = mod_par,
  var_par = cur_var_par,
  species_par = species_par,
  switch = switch,
  sim_id = sim_id,
  seed = master_seed,
  record_steps = steps_to_record
)

# Extract metadata from the last recorded state, to use for filenaming and logging.
meta <- results[[length(results)]]$meta

scenario <- scenario_key(meta = meta)
replicate_num <- rep_number(meta = meta)
filename <- sim_filename(sim_id, scenario, replicate_num)

cat("\nSimulation completed.\n\n")
if (export_states) {
  # We can now save all recorded steps into a single RDS file
  state_file <- paste0(state_dir, "/", filename, ".rds")
  saveRDS(results, file = state_file)
  cat(
    "Full states were recorded for time steps [",
    paste(names(results), collapse = ", "),
    "] and saved to:\n   ",
    state_file,
    "\n",
    sep = ""
  )
} else {
  cat("No states were exported. To export states, set export_states = TRUE.\n\n")
}


# Sampling ---------------------------------------------------------------

# Assess the recorded steps and prepare to sample the results.
recorded_steps <- names(results)
sampled_files <- character()

  # For each selected sampling method we will create a single CSV file 
  # with the sampled data from all recorded steps.
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

if (length(sampling_methods > 0)) {
  cat(
    "\nData was sampled for all time steps using methods [", 
    paste(sampling_methods, collapse = ", "), 
    "] and saved to:\n   ", 
    paste(sampled_files, collapse = "\n   "), 
    "\n", sep = ""
  )
} else {
  cat("No data was sampled. To sample data, specify at least one sampling method in sampling_methods.\n\n")
}

# Logging ----------------------------------------------------------------


# Finally, we can log the simulation details into the simulations log.
# The helper function automatically writes the log entry to the specified log file.
# Since every simulation has a unique sim_id, we shouldn't have concurrency issues.
if (to_log) {
  final_log <- log_entry(
    meta = meta,
    job_id = job_id, 
    scenario_key = scenario, 
    replicate_num = replicate_num,
    project_version = "2.0",
    status = "complete",
    state_file = state_file,
    sampled_files = sampled_files,
    log_file = log_file
  )
} else {
  cat("Simulation was not logged. To log the simulation, set to_log = TRUE.\n\n")
}
