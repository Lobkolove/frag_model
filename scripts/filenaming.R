# Script to fix all filenames of the output files to be consistent with the upcoming logging criteria.
# While fixing the names, we will also populate the log file with the relevant information for each simulation.
library(data.table)
library(here)
library(fs)
library(tools)
library(stringr)

old_state_dir <- "output/old_states"
old_sampled_dir <- "output/old_samples"

new_state_dir <- "output/model_states"
new_sampled_dir <- "output/sampled_data"

log_file <- "output/simulations_log.csv"

# Get all files in the state directory
state_files <- list.files(old_state_dir, full.names = TRUE)

# We will iterate over each state file, assess the relevant parameters and rename 
# all respective files (state + sampled) to match the new naming convention.

for (state_file in state_files) {

  # Extract old filename without extension
  old_filename <- file_path_sans_ext(basename(state_file))

  # Read the state file to extract parameters
  state_data <- readRDS(state_file)
  recorded_steps <- names(state_data)

  new_states <- vector("list", length(recorded_steps))

  for (step in recorded_steps) {
    
    single_state <- state_data[[step]]
    
    # Extract metadata 
    meta <- list(
      sim_id = single_state$sim_id,
      master_seed = single_state$master_seed,
      grid_size = single_state$grid_size,
      ac_amount = single_state$ac_amount,
      fragmentation = single_state$fragmentation,
      habitat = single_state$habitat,
      niche_breadth = 0.1,
      edge = 1, 
      dispersal = 1,
      dispersal_dist = 2
    )

    # Append reformated state to the new object
    new_states[[step]] <- list(
      meta = meta,
      grid = single_state$grid,
      agents = single_state$agents,
      ss_abund = single_state$ss_abund
    )
  }

  # Assess the replicate number (based on parameters) from the log file if it exists, otherwise assign 1
  if (file.exists(log_file)) {
    log <- fread(log_file)
    reps <- log[ac_amount == meta$ac_amount & fragmentation == meta$fragmentation & 
                  habitat == meta$habitat & niche_breadth == meta$niche_breadth & edge == meta$edge & 
                  dispersal == meta$dispersal & dispersal_dist == meta$dispersal_dist]
    replicate_num <- if (nrow(reps) > 0) max(reps$replicate_num, na.rm = TRUE) + 1 else 1
  } else {
    replicate_num <- 1
  }

  # Construct the new filename
  scenario <- paste0("ac", meta$ac_amount, "_frag", meta$fragmentation, "_edge", meta$edge, "_disp", meta$dispersal_dist)
  new_filename <- paste0(
    meta$sim_id,
    "_",
    scenario,
    "_r", sprintf("%03d", replicate_num)
  )

  # Save the new state file with the new name
  new_state_file <- file.path(new_state_dir, paste0(new_filename, ".rds"))
  saveRDS(new_states, new_state_file)

  # Now we will also rename the sampled files for this simulation
  sampling_methods <- list(all = "all", cb = "chessboard", rand = "random")
  new_sampled_files <- character()
      
  for (method in names(sampling_methods)) {
    
    old_sampled_file <- file.path(old_sampled_dir, paste0(old_filename, "_samp_", sampling_methods[[method]], ".csv"))
    table <- fread(old_sampled_file)

    new_sampled_file <- file.path(new_sampled_dir, paste0(new_filename, "_samp_", method, ".csv"))
    fwrite(table, new_sampled_file, na = "NA")
    new_sampled_files <- c(new_sampled_files, path_rel(new_sampled_file, start = here()))
    
  }

  # Finally, we will also update the log file with the relevant information for this simulation
  new_log_entry <- data.table(
    sim_id = meta$sim_id,
    job_id = "recovered", 
    scenario_key = scenario, 
    replicate_num = replicate_num,
    run_date = Sys.Date(), 
    project_version = "1.1",
    master_seed = meta$master_seed, 
    ac_amount = meta$ac_amount, 
    fragmentation = meta$fragmentation, 
    habitat = meta$habitat, 
    niche_breadth = meta$niche_breadth,
    dispersal = meta$dispersal, 
    dispersal_dist = meta$dispersal_dist, 
    edge = meta$edge,
    status = "complete",
    state_file = path_rel(new_state_file, start = here()),
    sampled_files = paste(new_sampled_files, collapse = "; ")
  )
  fwrite(new_log_entry, log_file, append = file.exists(log_file), logical01 = TRUE)

}