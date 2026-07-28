library(here)
library(data.table)
library(dplyr)
library(raster)
source(here("R", "sim_select.R"))

# Low frag
paths_lf <- sim_select(fragmentation == 0.2, file_type = "state")
state_low <- readRDS(paths_lf[1])
range(getValues(state_low$post_fragmentation$grid), na.rm = TRUE)
# 0.16 - 0.92

# min/max doesn't actually tell us how much of the whole range is actually covered,
# since the range is not necessarily contiguous
values_lf <- na.omit(round(getValues(state_low$post_fragmentation$grid), digits = 2))
hist(values_lf)
sort(unique(values_lf))
length(unique(values_lf)) / 100

# High frag
paths_hf <- sim_select(fragmentation == 0.8, file_type = "state")
state_high <- readRDS(paths_hf[1])
range(getValues(state_high$post_fragmentation$grid), na.rm = TRUE)
# 0.00 - 0.96

values_hf <- na.omit(round(getValues(state_high$post_fragmentation$grid), digits = 2))
hist(values_hf)
sort(unique(values_hf))
length(unique(values_hf)) / 100
