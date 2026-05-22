# Script to run an array job of simulations with different parameter combinations (cluster)
library(here)
library(data.table)
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
if (length(args) < 2) {
  stop("Usage: Rscript sim_pipeline.R <task_id> <job_id>")
}

task_id <- as.numeric(args[1])
job_id <- as.numeric(args[2])

# Paths
output_dir <- here("output")
log_file <- file.path(output_dir, "simulations_log.csv")
state_dir <- file.path(output_dir, "model_states")
sampled_dir <- file.path(output_dir, "sampled_data")

# Unique sim_id generation using block allocation to avoid concurrency issues
last_sim_id <- 0
existing_log <- NULL
if (file.exists(log_file)) {
  existing_log <- fread(log_file)
  last_sim_id <- max(as.numeric(existing_log$sim_id), 0, na.rm = TRUE)
}

# Block allocation: task_id gives offset within block
sim_id <- sprintf("%04d", last_sim_id + task_id)

cat("Block start =", last_sim_id, "→ sim_id =", sim_id, "task =", task_id, "\n")

# Parameter table
var_par_df <- tidyr::expand_grid(
  ac = 0.7,
  frag = c(0.2, 0.5, 0.8),
  hab = 0.15,
  nb = 0.1,
  disp = 1,
  disp_dist = c(1, 5, 10),
  edge = 1
)

# Determine which parameter combination to run based on task_id
n_scenarios <- nrow(var_par_df)
idx <- ((task_id - 1) %% n_scenarios) + 1
var_par <- var_par_df[idx, ]

# Compute scenario  key for logging
scenario <- paste0("ac", var_par$ac, "_frag", var_par$frag, "_edge", var_par$edge, "_disp", var_par$disp_dist)

# True per-parameter replicate numbering
replicate_num <- 1
if (!is.null(existing_log)) {
  reps <- existing_log[ac_amount == var_par$ac & fragmentation == var_par$frag & 
                        habitat == var_par$hab & niche_breadth == var_par$nb & edge_effect == var_par$edge & 
                        dispersal == var_par$disp & dispersal_dist == var_par$disp_dist]
  if (nrow(reps) > 0) replicate_num <- max(reps$replicate_num) + 1
}

cat("sim_id =", sim_id, "| task =", task_id, "| scenario =", scenario, 
    "| rep =", replicate_num, "| ac =", var_par$ac, "| frag =", var_par$frag, "\n")

# Run simulation
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

# Save model states
filename <- paste0(sim_id, "_", scenario, "_r", sprintf("%03d", replicate_num))
state_file <- file.path(state_dir, paste0(filename, ".rds"))
saveRDS(results, state_file)


# Sampling
recorded_steps <- names(results)
sampled_files <- character()
sampling_methods <- list(all = "all", cb = "checkerboard", rand = "random")

for (method in names(sampling_methods)) {
  sampled <- list()
  for (step in recorded_steps) {
    sampled[[step]] <- sample_cells(results[[step]], 
                                    method = sampling_methods[[method]], 
                                    n_samples = 30, 
                                    format = "long")
  }
  sampled_df <- rbindlist(sampled)
  
  samp_file <- file.path(sampled_dir, paste0(filename, "_samp_", method, ".csv")) 

  fwrite(sampled_df, samp_file, na = "NA")
  sampled_files <- c(sampled_files, samp_file)
}

# Single atomic log write to avoid concurrency issues
final_log <- data.table(
  sim_id = sim_id, 
  job_id = job_id, 
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
  edge_effect = var_par$edge,
  status = "complete",
  state_file = path_rel(state_file, start = here()),
  sampled_files = paste(sampled_files, collapse = "; ")
)

fwrite(final_log, log_file, append = file.exists(log_file), logical01 = TRUE)

cat("✓", sim_id, "\n\n")