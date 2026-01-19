require(raster)
require(data.table)
library(foreach)

source("Model/src/generate_grid.R")
source("Model/src/generate_agents.R")
source("Model/src/distribute_agents.R")
source("Model/src/birth.R")
source("Model/src/death.R")
source("Model/src/initialize.R")
source("Model/src/animate.R")
source("Model/src/cookie_cutting.R")
source("Model/src/immigration.R")
source("Model/src/landscape.R")
source("Model/src/disperse.R")
source("Model/src/run_model_step.R")
source("Model/src/clean_run.R")
source("Model/parameters.R")
source("R/sample_cells.R")


# Run dynamic model and get full state right after fragmentation
sim_id <- 1
states <- clean_run(
  mod_par,
  var_par,
  switch,
  sim_id,
  record_steps = "post_frag_start",
  seed
) 

post_frag_state <- states$post_frag_start

# Export as RDS
saveRDS(post_frag_state, "data-raw/post_frag_state_test1.rds")


# Sampling ----------------------------------------------------------------

# 30 random habitat cells
sample_ran <- sample_cells(post_frag_state, method = "random", n_samples = 30)

# All habitat cells
sample_all <- sample_cells(post_frag_state, method = "all")

# Every other cell (chessboard)
sample_cb <- sample_cells(post_frag_state, method = "chessboard")

# Export to CSVs
write.csv(sample_ran, "data-raw/post_frag_sample_ran.csv")
write.csv(sample_all, "data-raw/post_frag_sample_full.csv")
write.csv(sample_cb, "data-raw/post_frag_sample_cb.csv")




