# Reduced version of the filenaming script which only adds the missing information to the log file
library(data.table)
library(here)
library(tools)

source("Model/src/registry.R")

states_dir <- "output/model_states"
sampled_dir <- "output/sampled_data"
log_file <- "output/simulations_log.csv"

# Get all state files
state_files <- list.files(states_dir, full.names = TRUE)

# Loop through each state file and extract metadata to update the log
for (state_file in state_files) {

  # Extract base filename without extension
  base_filename <- file_path_sans_ext(basename(state_file))

  # Read the state file to extract parameters
  state_data <- readRDS(state_file)
  recorded_steps <- names(state_data)

  # We use post_fragmentation step to extract the parameters, 
  # since it is always recorded and contains all relevant metadata
  single_state <- state_data[["post_fragmentation"]]
  
  # Currently, all saved states should have a meta object, but some old simulations didn't have it. 
  # So just to be sure, we will check if it exists, and if not, we will build it from the state parameters.
  if (is.null(single_state$meta)) {
    # Extract metadata
    meta <- list(
      sim_id = single_state$sim_id,
      master_seed = single_state$master_seed,
      grid_size = single_state$grid_size,
      ac_amount = single_state$ac_amount,
      habitat = single_state$habitat,
      fragmentation = single_state$fragmentation,
      # All other relevant parameters were fixed in the old simulations, so we can just assign them here
      niche_breadth = 0.1,
      edge_effect = 1,
      dispersal = 1,
      dispersal_dist = 2
    )
  } else {
    meta <- single_state$meta
    # In some cases edge_effect was called edge
    if (is.null(meta$edge_effect)) meta$edge_effect <- meta$edge
  }

  # We now assess if there are sampled files for this simulation (sharing the same base filename)
  sampled_files <- list.files(
    path = sampled_dir,
    pattern = base_filename,
    full.names = TRUE
  )

  # Add log entry for this simulation
  log_entry(
    meta = meta,
    job_id = "recovered",
    scenario_key = scenario_key(meta = meta),
    replicate_num = rep_number(meta = meta),
    project_version = "1.2",
    status = "complete",
    state_file = state_file,
    sampled_files = sampled_files,
    overwrite = FALSE # Should entries be overwritten if already present?
  )
}
