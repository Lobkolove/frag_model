library(data.table)
library(here)
library(digest)
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

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript sim_pipeline.R <task_id> [job_id]")
}

task_id <- as.numeric(args[1])
job_id <- if (length(args) > 1) args[2] else "local"

# Paths
output_dir <- here("output")
log_file <- file.path(output_dir, "simulations_log.csv")
state_dir <- file.path(output_dir, "model_states")
sampled_dir <- file.path(output_dir, "sampled_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(state_dir, showWarnings = FALSE)
dir.create(sampled_dir, showWarnings = FALSE)

# Minimal sim_id (4-digit, auto-increment)
next_sim_id <- 1
if (file.exists(log_file)) {
  for (i in 1:5) {  # Retry with backoff
    tryCatch({
      existing_log <- fread(log_file)
      next_sim_id <- max(as.numeric(existing_log$sim_id), 1, na.rm = TRUE) + 1
      break
    }, error = function(e) {
      Sys.sleep(0.1 * i)
      if (i == 5) stop("Cannot read log: ", e$message)
    })
  }
}
sim_id <- sprintf("%04d", next_sim_id)

# Parameter table
var_par_df <- tidyr::expand_grid(
  ac = seq(0.1, 0.9, 0.2), frag = c(0.2, 0.5, 0.8),
  hab = 0.15, nb = 0.1, disp = 1, disp_dist = 2, edge = 1
)
var_par <- var_par_df[task_id, ]

# Compute scenario  key for logging
scenario_key <- paste0("ac", sprintf("%.1f", var_par$ac), "_frag", var_par$frag, "_hab", sprintf("%.2f", var_par$hab))

# True per-parameter replicate numbering
replicate_num <- 1
if (file.exists(log_file)) {
  log <- fread(log_file)
  scenario_reps <- log[scenario_key == scenario_key]
  if (nrow(scenario_reps) > 0) {
    replicate_num <- max(scenario_reps$replicate_num, 0) + 1
  }
}

cat("sim_id =", sim_id, "- task =", task_id, "- scenario =", scenario_key, 
    "- rep =", replicate_num, "- ac =", var_par$ac, "frag =", var_par$frag, "\n")

# Run simulation
seed_used <- master_seed + task_id
results <- clean_run(
  mod_par = mod_par,
  var_par = var_par,
  switch = switch,
  sim_id = sim_id,
  seed = seed_used,
  record_steps = c(
    # "start",
    # "pre_fragmentation",
    "post_fragmentation",
    "final"
  )
)

# Save model states
recorded_steps <- names(results)
for (step in recorded_steps) {
  step_file <- file.path(state_dir, paste0(sim_id, "_", scenario_key, 
                                          "_r", sprintf("%03d", replicate_num), 
                                          "_", step, ".rds"))
  saveRDS(results[[step]], step_file, compress = "xz")
}

# Sampling
sampled_files <- character()
sampling_methods <- c("all", "chessboard", "random")

for (method in sampling_methods) {
  sampled <- list()
  for (step in recorded_steps) {
    sampled[[step]] <- sample_cells(results[[step]], 
                                    method = method, 
                                    n_samples = 30, 
                                    format = "long")
  }
  sampled_df <- rbindlist(sampled)
  sampled_df[, `:=`(replicate_num = replicate_num)]
  
  samp_file <- file.path(sampled_dir, paste0(sim_id, "_", scenario_key, 
                                            "_r", sprintf("%03d", replicate_num), 
                                            "_samp_", method, ".csv"))
  fwrite(sampled_df, samp_file, na = "NA")
  sampled_files <- c(sampled_files, samp_file)
}

# Single atomic log write to avoid concurrency issues
final_log <- data.table(
  sim_id, task_id, job_id, scenario_key, replicate_num,
  run_date = Sys.Date(), project_version = "frag_v1",
  ac = var_par$ac, frag = var_par$frag, hab = var_par$hab, nb = var_par$nb,
  disp = var_par$disp, disp_dist = var_par$disp_dist, edge = var_par$edge,
  seed = seed_used, status = "complete",
  state_files = paste(state_files, collapse = "; "), 
  sampled_files = paste(sampled_files, collapse = "; ")
)

temp_file <- tempfile(pattern = "log_", fileext = ".csv")
if (file.exists(log_file)) {
  existing <- fread(log_file)
  fwrite(rbind(existing, final_log, fill = TRUE), temp_file)
} else {
  fwrite(final_log, temp_file)
}
file.rename(temp_file, log_file)

cat("✓", sim_id, "\n\n")