# Script to fix all filenames of the output files to be consistent with the upcoming logging criteria.
# While fixing the names, we will also populate the log file with the relevant information for each simulation.
library(data.table)
library(here)
library(fs)
library(tools)
library(stringr)

source("Model/src/registry.R")
source("Model/src/record_state.R")

old_state_dir <- "output/old_states"
old_sampled_dir <- "output/old_samples"

new_state_dir <- "output/model_states"
new_sampled_dir <- "output/sampled_data"

# Create new directories if they don't exist
dir_create(new_state_dir)
dir_create(new_sampled_dir)

log_file <- "output/simulations_log.csv"

# Get all files in the state directory
state_files <- list.files(old_state_dir, full.names = TRUE)

# We will iterate over each state file, assess the relevant parameters and rename
# all respective files (state + sampled) to match the new naming convention.

for (state_file in state_files) {
  # Extract old filename without extension
  old_filename <- tools::file_path_sans_ext(basename(state_file))

  # Read the state file to extract parameters
  state_data <- readRDS(state_file)
  recorded_steps <- names(state_data)

  # Initialize list vector with length and names of recorded steps
  new_states <- vector("list", length(recorded_steps))
  names(new_states) <- recorded_steps

  for (step_label in recorded_steps) {
    single_state <- state_data[[step_label]]

    # Check if state has a meta object
    if (is.null(single_state$meta)) {
      # Construct metadata
      meta <- build_meta(
        sim_id = single_state$sim_id,
        run_date = as.Date(file.info(state_file)$ctime),
        master_seed = single_state$master_seed,
        grid_size = single_state$grid_size,
        ac_amount = single_state$ac_amount,
        habitat = single_state$habitat,
        fragmentation = single_state$fragmentation,
        niche_breadth = 0.1,
        edge_effect = 1,
        dispersal_type = "short_long",
        dispersal_kernel = "exponential",
        dispersal = 1,
        dispersal_dist = 2
      )
    } else {
      # Extract and standardise metadata
      meta <- build_meta(single_state$meta)
    }

    # Append reformatted state to the new object
    new_states[[step_label]] <- record_state(
      model_state = single_state,
      meta = meta,
      step_label = step_label
    )
  }

  # Assess the replicate number (based on parameters) from the log file if it exists, otherwise assign 1
  replicate_num <- rep_number(meta)

  # Construct the new filename
  scenario <- scenario_key(meta)
  new_filename <- sim_filename(
    sim_id = meta$sim_id,
    scenario_key = scenario,
    replicate_num = replicate_num
  )

  # Save the new state file with the new name
  new_state_file <- file.path(new_state_dir, paste0(new_filename, ".rds"))
  saveRDS(new_states, new_state_file)

  # Now we will also rename the sampled files for this simulation
  sampling_methods <- list(all = "all", cb = "chessboard", rand = "random")
  new_sampled_files <- character()

  for (method in names(sampling_methods)) {

    old_sampled_file <- list.files(
      old_sampled_dir,
      pattern = paste0(old_filename, ".*", method, "\\.csv$"),
      full.names = TRUE
    )
    
    if (length(old_sampled_file) == 0) {
      old_sampled_file <- list.files(
        old_sampled_dir,
        pattern = paste0(old_filename, ".*", sampling_methods[[method]], "\\.csv$"),
        full.names = TRUE
      )
    }

    table <- fread(old_sampled_file)

    # Add missing columns and reorder columns to match the new log format
    table <- table %>%
      dplyr::mutate(
        habitat = meta$habitat,
        niche_breadth = meta$niche_breadth,
        edge_effect = meta$edge_effect,
        dispersal_type = meta$dispersal_type,
        dispersal_ratio = meta$dispersal_ratio,
        dispersal_dist = meta$dispersal_dist
      ) %>%
      dplyr::select(
        sim_id,
        master_seed,
        grid_size,
        fragmentation,
        ac_amount,
        habitat,
        niche_breadth,
        edge_effect,
        dispersal_type,
        dispersal_ratio,
        dispersal_dist,
        step,
        step_label,
        samp_method,
        sample_id,
        cell_id,
        x_loc,
        y_loc,
        patch_id,
        patch_size,
        species_id,
        n
      )

    new_sampled_file <- file.path(
      new_sampled_dir,
      paste0(new_filename, "_", method, ".csv")
    )
    fwrite(table, new_sampled_file, na = "NA")
    new_sampled_files <- c(
      new_sampled_files,
      path_rel(new_sampled_file, start = here())
    )
  }

  # Finally, we will also update the log file with the relevant information for this simulation
  log_entry(
    meta = meta,
    job_id = "recovered",
    scenario_key = scenario,
    replicate_num = replicate_num,
    project_version = "1.1",
    status = "complete",
    state_file = new_state_file,
    sampled_files = new_sampled_files,
    overwrite = TRUE
  )
}
  