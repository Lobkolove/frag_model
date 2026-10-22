################################################################################
#
# Some examples of how to run the new model (Selina)
# 1. Just one model run
# 2. Multiple runs in a loop
# 3. Multiple loops in parallel
################################################################################

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


# Setup -------------------------------------------------------------------

# set the numbers of repetitions
task_id <- 1


# Single model run --------------------------------------------------------
# one model repetiton for the first parameter set
result <- GeDo_run(
  mod_par = mod_par,
  var_par = var_par[1, ],
  switch = switch,
  file_name = "test_Selina",
  task_id = 1,
  sim_id = 1
)

# save the results
data.table::fwrite(
  x = result[["output_all"]],
  file = paste0(
    "Outputs/heatmap",
    test_name, "_rep_", task_id, "_output_general.csv"
  ),
  sep = ","
)

data.table::fwrite(
  x = result[["output_sample"]],
  file = paste0(
    "Outputs/heatmap",
    test_name, "_rep_", task_id, "_output_sample.csv"
  ),
  sep = ","
)


# Run in a loop (not parallel) --------------------------------------------
result_seq <- foreach(
  i = 1:5
) %do% {
  GeDo_run(
    mod_par = mod_par,
    var_par = var_par[i, ],
    switch = switch,
    file_name = "test2",
    task_id = 1,
    sim_id = i
  )
}

# transpose the list, so all output_all and all output_sample are togehter
result_par <- purrr::transpose(result_seq)

# rbind all results
output_all <- data.table::rbindlist(
  result_par[["output_all"]], 
  idcol = "sim_id"
)
output_sample <- data.table::rbindlist(
  result_par[["output_sample"]], 
  idcol = "sim_id"
)

# write them to disk
data.table::fwrite(
  x = output_all,
  file = paste0(
    "Outputs/heatmap",
    test_name, "_rep_", task_id, "_output_general.csv"
  ),
  sep = ","
)

data.table::fwrite(
  x = output_sample,
  file = paste0(
    "Outputs/heatmap",
    test_name, "_rep_", task_id, "_output_sample.csv"
  ),
  sep = ","
)

# Run in a loop (parallel) ------------------------------------------------
# how many cores to use (here 1 less than is available)
workers <- parallel::detectCores() - 1
my_cluster <- parallel::makeCluster(
  workers,
  type = "PSOCK" # PSOCK works on both Win and UNIX
)
doParallel::registerDoParallel(cl = my_cluster)

result_par <- foreach(
  i = 1:5,
  .packages = c("data.table", "raster")
) %dopar% {
  GeDo_run(
    mod_par = mod_par,
    var_par = var_par[i, ],
    switch = switch,
    file_name = "test2",
    task_id = 1,
    sim_id = i
  )
}

# Combine and write results -----------------------------------------------
# transpose the list, so all output_all and all output_sample are togehter
result_par <- purrr::transpose(result_par)

# rbind all results
output_all <- data.table::rbindlist(result_par[["output_all"]], idcol = "sim_id")
output_sample <- data.table::rbindlist(result_par[["output_sample"]], idcol = "sim_id")

# write them to disk
data.table::fwrite(
  x = output_all,
  file = paste0(
    "Outputs/heatmap",
    test_name, "_rep_", task_id, "_output_general.csv"
  ),
  sep = ","
)

data.table::fwrite(
  x = output_sample,
  file = paste0(
    "Outputs/heatmap",
    test_name, "_rep_", task_id, "_output_sample.csv"
  ),
  sep = ","
)
