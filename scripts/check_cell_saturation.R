# How strongly did cell carrying capacity constrain the communities at the
# pre-fragmentation, post-fragmentation and final time steps?
#
# Capacity never enters mortality: death() draws independently per individual at
# death_rate, and k_inter/k_intra are only checked in birth(), immigration() and
# distribute_agents(). What capacity does is block recruitment.
#
# States are recorded after death (birth -> death -> immigration -> record in
# run_model_step()), so a cell that filled up to k_inter during the birth loop is
# recorded only after roughly death_rate of its occupants were killed. A cell
# sitting at capacity therefore shows up at about k_inter * (1 - death_rate)
# individuals, i.e. 37.5 rather than 50, with Binomial survival noise around it.
# That reference value is what cell abundances are scaled against below, so a
# saturation index of 1 means "at capacity" and a plain count above some fixed
# threshold is not needed.
#
# For the moment this check covers the core series only.

library(here)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)

source("R/sim_select.R")

# Plot settings
theme_set(theme_bw(base_size = 14))
pal_frag <- c("#4688ad", "#e9b14a", "#ac5384")

# Model parameters the reference value is built from (see Model/parameters.R)
k_inter <- 50
death_rate <- 0.25

# Expected number of individuals recorded in a cell that was at capacity
expected_at_cap <- k_inter * (1 - death_rate)

# Time steps to check
steps <- c("pre_fragmentation", "post_fragmentation", "final")

# Read in log file
log_file <- here("output/simulations_log.csv")
log <- fread(log_file)


# 1. Core series ---------------------------------------------------------

# Assess core simulation IDs
core_ids <- sim_ids(habitat == 0.15, ac_amount == 0.7, dispersal_dist == 2)

# Check that no simulations share the same master seed
core_seeds <- log %>%
  filter(sim_id %in% core_ids) %>%
  pull(master_seed)
sum(duplicated(core_seeds))
# No duplicated master seeds, ready for analysis!

# Select simulation and get path for full sample file
paths_core_full <- sim_select(sim_id %in% core_ids, sampled = "all")

# Read in data, keeping only the three steps of interest
data_core_full <- map(
  here(paths_core_full),
  ~ fread(.x) %>%
    dplyr::filter(step_label %in% steps)
)


  ## 1.1 Cell saturation ---------------------------------------------------

#' Measure how close cells sit to their carrying capacity
#'
#' Sums individuals across species within each cell and step, then scales the
#' mean cell abundance against the abundance expected of a cell that was at
#' capacity when births were drawn. Empty cells carry no information about
#' capacity, so the index is computed over occupied cells only and occupancy is
#' reported separately.
#'
#' @param data A sampled-data table of a single simulation.
#' @param expected_at_cap Abundance expected of a cell recorded at capacity,
#'   i.e. `k_inter * (1 - death_rate)`.
#' @param steps Time steps to report on, in the order they should be shown.
#'
#' @return One row per simulation and step, with the saturation index over
#'   occupied cells, the same index over all cells, occupancy, and the largest
#'   cell abundance observed.
cell_saturation <- function(data, expected_at_cap = 37.5, steps = c("pre_fragmentation", "post_fragmentation", "final")) {

  data %>%
    dplyr::filter(step_label %in% steps) %>%
    dplyr::group_by(sim_id, fragmentation, step_label, cell_id) %>%
    dplyr::summarise(n_ind = sum(n, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(sim_id, fragmentation, step_label) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_occupied = sum(n_ind > 0),
      prop_occupied = n_occupied / n_cells,
      mean_n_occ = mean(n_ind[n_ind > 0]),
      sat_index = mean_n_occ / expected_at_cap,
      sat_index_all = mean(n_ind) / expected_at_cap,
      max_n_ind = max(n_ind),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      step_label = factor(step_label, levels = steps)
    )
}

# Measure saturation for each dataset in the core series
sat_core_full <- map(
  data_core_full,
  ~ cell_saturation(.x, expected_at_cap = expected_at_cap, steps = steps)
) %>%
  bind_rows() %>%
  mutate(
    fragmentation = factor(
      fragmentation,
      levels = c(0.2, 0.5, 0.8),
      labels = c("Low", "Medium", "High")
    ),
    step_label = factor(
      step_label,
      levels = steps,
      labels = c("Pre-fragmentation", "Post-fragmentation", "End of simulation")
    )
  )

# Per-simulation saturation
print(sat_core_full, n = nrow(sat_core_full))

# Note: fragmentation is only defined once the landscape is fragmented, so the
# pre-fragmentation step carries NA and is summarised across all replicates.
sat_core_summary <- sat_core_full %>%
  group_by(step_label, fragmentation) %>%
  summarise(
    n_sims = n(),
    mean_sat = mean(sat_index),
    q_low = quantile(sat_index, probs = 0.025),
    q_high = quantile(sat_index, probs = 0.975),
    mean_sat_all = mean(sat_index_all),
    mean_occupied = mean(prop_occupied),
    max_n_ind = max(max_n_ind),
    .groups = "drop"
  )
sat_core_summary


  ## 1.2 Plots -------------------------------------------------------------

    ### 1.2.1 Saturation by time step -----

# Fragmentation levels pooled, one box per time step.
# The dashed line marks a cell sitting exactly at carrying capacity.
gg_sat_core_step <- ggplot(sat_core_full, aes(step_label, sat_index)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_boxplot(fill = "grey85", width = 0.5, outlier.alpha = 0.5) +
  labs(
    x = "Time step",
    y = "Saturation index (mean occupancy / capacity)"
  )
gg_sat_core_step

    ### 1.2.2 Saturation by fragmentation -----

# Fragmentation is only defined once the landscape is fragmented, so the
# pre-fragmentation step is dropped here.
gg_sat_core_frag <- sat_core_full |>
  dplyr::filter(step_label != "Pre-fragmentation") |>
  ggplot(aes(fragmentation, sat_index, fill = fragmentation)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.alpha = 0.5) +
  facet_wrap(~ step_label) +
  labs(
    x = "Level of fragmentation",
    y = "Saturation index (mean occupancy / capacity)",
    fill = "Level of\nfragmentation"
  ) +
  scale_fill_manual(values = pal_frag) +
  theme(legend.position = "none")
gg_sat_core_frag
