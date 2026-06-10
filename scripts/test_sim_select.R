library(data.table)
library(dplyr)
library(stringr)
source("R/sim_select.R")

sim_select

# All state files from simulations with varying fragmentation and all other parameters as static
sim_select(file_type = "state")

# All sampled files from simulations with varying fragmentation AND dispersal distance
sim_select(
  vars = c("fragmentation", "dispersal_dist"),
  file_type = "sampled",
  sampled = "random"
)
# Note that by default the function selects the parameter combination with the biggest number of simulations.
# If you want to select a specific parameter combination, you can set `mode` to "user_defined" 
sim_select(
  vars = c("fragmentation", "dispersal_dist"),
  file_type = "sampled",
  sampled = "random",
  mode = "user_defined"
)

# The function can also take specific conditions (e.g. fragmentation = 0.5 and dispersal_dist = 10)
sim_select(
  fragmentation = 0.5,
  dispersal_dist = 10,
  file_type = "sampled",
  sampled = "random"
)
