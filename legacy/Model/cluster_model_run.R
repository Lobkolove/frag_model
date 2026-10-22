################################################################################
#
# Preparation - Loading required packages and functions
#
################################################################################
# 
library(raster)
library(data.table)
library(tibble)
library(foreach)

source("src/generate_grid.R")
source("src/generate_agents.R")
source("src/distribute_agents.R")
source("src/multiple_runs.R")
source("src/birth.R")
source("src/death.R")
source("src/initialize.R")
source("src/run_dynamic_model.R")
source("src/multiple_run_dynamic_model.R")
source("src/animate.R")
source("src/GeDo_run.R")
source("src/cookie_cutting.R")
source("src/immigration.R")
source("src/landscape.R")
source("src/disperse.R")
source("parameters.R")


################################################################################
#
# Model Run
#
################################################################################

# set the numbers of repetitions

args <- commandArgs(trailingOnly = TRUE)
task_id <- args[1]
workers <- as.integer(args[2])
exp_name <- "heatmap_ac_"
exp_num <- "29" # don't forget to change!

paste0("I am running with ", workers, " workers")

paste0("And ", parallel::detectCores(), " cores")

# prepare parallel session ------------------------------------------------

my.cluster <- parallel::makeCluster(
  workers,
  port = 11000 + as.integer(task_id),
  outfile = "",
  type = "FORK" # Change to "PSOCK" if running locally!
)

# check cluster definition (optional)
print("Cluster definition ------------")
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)

# check if it is registered (optional)
print("Check registered cluster ------------")
foreach::getDoParRegistered()

# how many workers are available? (optional)
print("Check workers ------------")
foreach::getDoParWorkers()

result_par <- foreach(
  i = 1:nrow(var_par),
  .packages = c("data.table", "raster")
) %dopar% {
  GeDo_run(
    mod_par = mod_par,
    var_par = var_par[i, ],
    switch = switch,
    file_name = exp_name,
    task_id = 1,
    sim_id = i
  )
}

parallel::stopCluster(cl = my.cluster)

# Combine and write results -----------------------------------------------

# create static output files
static_par <- rbind(as.data.frame(t(switch)), as.data.frame(t(mod_par)))
write.csv(static_par, file = paste0("Outputs/", exp_name, exp_num, "_static_parameters.csv"))

var_par_write <- var_par
var_par_write <- rowid_to_column(var_par_write, "sim_ID")

write.csv(var_par_write, file = paste0("Outputs/", exp_name, exp_num, "_varaying_parameters.csv"), row.names = FALSE)

# transpose the list, so all output_all and all output_sample are together
result_par <- purrr::transpose(result_par)

# rbind all results
output_all <- data.table::rbindlist(result_par[["output_all"]])
output_sample <- data.table::rbindlist(result_par[["output_sample"]])

# write them to disk
data.table::fwrite(
  x = output_all,
  file = paste0(
    "Outputs/",
    exp_name, exp_num, "_rep_", task_id, "_output_general.csv"
  ),
  sep = ","
)

data.table::fwrite(
  x = output_sample,
  file = paste0(
    "Outputs/",
    exp_name, exp_num, "_rep_", task_id, "_output_sample.csv"
  ),
  sep = ","
)
