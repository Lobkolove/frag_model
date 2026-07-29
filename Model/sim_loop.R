# Script for a simulation loop (local)
# Can be used for a single run by keeping var_par as a single row data frame (see parameters.R).
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

# Paths
log_file <- "output/simulations_log.csv"
state_dir <- "output/model_states"
sampled_dir <- "output/sampled_data"

# Notes (any additional information to be printed to the console output)
notes <- "k_intra = 5"


# Settings ---------------------------------------------------------------

# Select the time steps to be recorded and sampled.
steps_to_record <- c(
  # "start",
  # "pre_fragmentation",
  "post_fragmentation",
  "final"
)
# Note that if steps_to_record is empty, simulations will still run but
# no states will be recorded and no sampling can be performed.

# Should full states be exported as RDS files? 
# If yes, set export_states = TRUE.
export_states <- FALSE

# Should the results be sampled? 
# If yes, specify the sampling method(s) by uncommenting the desired method(s) below.
# If no, assign sampling_methods as empty list().
sampling_methods <- list(
  # all = "all", 
  # cb = "checkerboard", 
  # rand = "random"
)

# Should the simulations be logged? 
# If yes, set to_log = TRUE. 
to_log <- FALSE
# All simulations relevant for analyses should be logged!
# Set to_log = FALSE only for testing and debugging purposes!

# Loop -------------------------------------------------------------

# Assess the number of simulations to run based on the parameter grid (var_par)
n_sims <- nrow(var_par)
results_list <- vector(mode = "list", length = n_sims)

# Simulations loop
for (i in seq_len(n_sims)) {
  
  ## Simulation -------------------------------------------------------------
  
  # Define a unique identifier for this simulation run (sim_id)
  if (to_log) {
    # Last used sim_id is determined from the log file, and the next sim_id is assigned.
    sim_id <- unique_sim_id(log_file)
  } else {
    # If not logging, we can just use the loop index as sim_id
    sim_id <- i
  }

  # A single master seed is set in parameters.R.
  # If we want to have different seeds for each simulation, 
  # we need to increment the master_seed for each simulation.
  master_seed <- master_seed + i - 1

  # Extract parameters for this simulation
  cur_var_par <- var_par[i, ]
  
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

  # Store the results for this simulation
  results_list[[i]] <- results

  # Extract metadata from the last recorded state, to use for filenaming and logging.
  meta <- results[[length(results)]]$meta
  
  # Similarly to what we did with sim_id, we can also assign a replicate number for this scenario  
  scenario <- scenario_key(meta = meta)
  replicate_num <- rep_number(meta = meta)
  filename <- sim_filename(sim_id, scenario, replicate_num)
  
  cat("\nSimulation completed.\n\n ")
  if (export_states) {
    # We can now save all recorded steps into a single RDS file
    state_file <- paste0(state_dir, "/", filename, ".rds")
    saveRDS(results, state_file)
    cat(
      "\nFull states were recorded for time steps [", 
      paste(names(results), collapse = ", "), 
      "] and saved to:\n   ", state_file, "\n", 
      sep = ""
    )
  } else {
  cat("No states were exported. To export states, set export_states = TRUE.\n\n")
  }
  
  ## Sampling ---------------------------------------------------------------
  
  # Assess the recorded steps and prepare to sample the results.
  recorded_steps <- names(results)
  sampled_files <- character()
  
  # For each selected sampling method we will create a single CSV file 
  # with the sampled data from all recorded steps.
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
  
  ## Logging ----------------------------------------------------------------
  
  # Finally, we can log the simulation details into the simulations log.
  # This is important for keeping track of what simulations we have run, the parameters used, 
  # and where the outputs are saved. 
  # It also helps to avoid overwriting existing simulations and to identify which replicates belong to which scenarios.
  if (to_log) {
    log_entry(
      meta = meta,
      job_id = "local", 
      scenario_key = scenario, 
      replicate_num = replicate_num,
      project_version = "2.0",
      status = "complete",
      state_file = state_file,
      sampled_files = sampled_files
    )
  } else {
    cat("Simulation was not logged. To log the simulation, set to_log = TRUE.\n\n")
  }
}
