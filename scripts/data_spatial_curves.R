require(raster)
require(data.table)
library(foreach)

source("Model/src/generate_grid.R")
source("Model/src/generate_agents.R")
source("Model/src/distribute_agents.R")
source("Model/src/multiple_runs.R")
source("Model/src/birth.R")
source("Model/src/death.R")
source("Model/src/initialize.R")
source("Model/src/run_dynamic_model.R")
source("Model/src/multiple_run_dynamic_model.R")
source("Model/src/animate.R")
source("Model/src/GeDo_run.R")
source("Model/src/cookie_cutting.R")
source("Model/src/immigration.R")
source("Model/src/landscape.R")
source("Model/src/disperse.R")
source("Model/parameters.R")


# Low fragmentation:
var_par$frag <- 0.2
result <- GeDo_run(
  mod_par = mod_par,
  var_par = var_par[1, ],
  switch = switch,
  file_name = "test_Ju",
  task_id = 1,
  sim_id = 1
)

data.table::fwrite(
  x = result[["output_sample"]],
  file = paste0("data-raw/test_Ju_low_rep_", task_id, "_output_sample.csv"),
  sep = ","
)

# High fragmentation:
var_par$frag <- 0.8
result <- GeDo_run(
  mod_par = mod_par,
  var_par = var_par[1, ],
  switch = switch,
  file_name = "test_Ju",
  task_id = 1,
  sim_id = 1
)

data.table::fwrite(
  x = result[["output_sample"]],
  file = paste0("data-raw/test_Ju_high_rep_", task_id, "_output_sample.csv"),
  sep = ","
)
