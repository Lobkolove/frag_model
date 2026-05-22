library(data.table)
library(tools)
source("R/sample_cells.R")
source("R/toroidal_clump.R")

old_sampled_dir <- "output/old_samples"
old_states_dir <- "output/old_states"

# List all state files starting with 0076 - 0084
state_files <- list.files(old_states_dir, pattern = "^007[6-9]|008[0-4]", full.names = TRUE)

sampling_methods <- c(all = "all", checkerboard = "chessboard", random = "random")

# Loop through each state file, sample cells, and save results
for (state_file in state_files) {

  filename <- file_path_sans_ext(basename(state_file))

  # Read in states
  states <- readRDS(state_file)

  recorded_steps <- names(states)

  for (method in names(sampling_methods)) {
    sampled <- list()
    for (step in recorded_steps) {
      sampled[[step]] <- sample_cells(states[[step]], 
                                      method = method, 
                                      n_samples = 30, 
                                      format = "long")
    }
    sampled_df <- rbindlist(sampled)

    samp_file <- file.path(old_sampled_dir, paste0(filename, "_samp_", sampling_methods[[method]], ".csv"))
    fwrite(sampled_df, samp_file, na = "NA")
  }

}
