# Mean alpha/gamma species richness ratios for the main simulation series
#
# Uses the same simulation sets and sample handling as scripts/main_analyses.R:
# full ("_all") samples, richness at the end of the simulation, and for the
# dispersal series the same exclusion of replicate 006 (datasets 46–54).
#
# Output: output/exploration/alpha_gamma_ratio_final.csv

library(here)
library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

source("R/sim_select.R")
source("R/compute_diversity.R")

steps <- c("post_fragmentation", "final")

log_file <- here("output/simulations_log.csv")
log <- fread(log_file)

read_series_data <- function(ids) {
  paths <- sim_select(sim_id %in% ids, sampled = "all")
  map(
    here(paths),
    ~ fread(.x) |>
      dplyr::filter(step_label %in% steps)
  )
}

#' Per-simulation alpha (mean sample richness) and gamma (landscape richness)
compute_series_diversity <- function(data_list, metadata_cols = NULL) {
  bind_rows(
    map(
      data_list,
      ~ compute_diversity(data = .x, metadata_cols = metadata_cols)
    )
  )
}

#' Collapse to scenario × fragmentation summaries used in the exploration CSV
summarize_alpha_gamma_ratios <- function(
  div_long,
  series,
  scenario_col = NULL,
  scenario_label = NULL
) {
  div_long <- div_long |>
    filter(index == "richness", step_label == "End of simulation")

  if (is.null(scenario_col)) {
    per_sim <- div_long |>
      mutate(scenario = scenario_label) |>
      select(sim_id, scenario, fragmentation, scale, value)
  } else {
    per_sim <- div_long |>
      mutate(scenario = paste0(scenario_col, " = ", .data[[scenario_col]])) |>
      select(sim_id, scenario, fragmentation, scale, value, all_of(scenario_col))
  }

  per_sim <- per_sim |>
    pivot_wider(names_from = scale, values_from = value) |>
    mutate(ratio = `Local scale` / `Landscape scale`)

  per_sim |>
    group_by(scenario, fragmentation) |>
    summarise(
      n_sims = n(),
      mean_alpha = mean(`Local scale`, na.rm = TRUE),
      mean_gamma = mean(`Landscape scale`, na.rm = TRUE),
      ratio_of_means = mean_alpha / mean_gamma,
      mean_of_ratios = mean(ratio, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(series = series) |>
    select(
      series,
      scenario,
      fragmentation,
      n_sims,
      mean_alpha,
      mean_gamma,
      ratio_of_means,
      mean_of_ratios
    )
}


# 1. Core series ---------------------------------------------------------

core_ids <- sim_ids(habitat == 0.15, ac_amount == 0.7, dispersal_dist == 2)
div_core <- compute_series_diversity(read_series_data(core_ids))
out_core <- summarize_alpha_gamma_ratios(
  div_core,
  series = "core",
  scenario_label = "core"
)


# 2. Habitat amount series -----------------------------------------------

habitat_ids <- sim_ids(habitat %in% c(0.1, 0.3, 0.5))
div_hab <- compute_series_diversity(
  read_series_data(habitat_ids),
  metadata_cols = "habitat"
)
out_hab <- summarize_alpha_gamma_ratios(
  div_hab,
  series = "habitat",
  scenario_col = "habitat"
)


# 3. Dispersal distance series -------------------------------------------

dispersal_ids <- sim_ids(dispersal_dist %in% c(1, 4, 8))
dispersal_data <- read_series_data(dispersal_ids)
# Same replicate excluded in main_analyses.R (8 cells / high fragmentation)
dispersal_data[46:54] <- NULL

div_disp <- compute_series_diversity(
  dispersal_data,
  metadata_cols = "disp_dist"
)
out_disp <- summarize_alpha_gamma_ratios(
  div_disp,
  series = "dispersal",
  scenario_col = "disp_dist"
)


# 4. Environmental autocorrelation series ------------------------------

ac_ids <- sim_ids(ac_amount %in% c(0, 0.5, 1), dispersal_type == "short_long")
ac_seeds <- log |>
  filter(sim_id %in% ac_ids) |>
  pull(master_seed)
ac_ids <- ac_ids[!duplicated(ac_seeds)]

div_ac <- compute_series_diversity(
  read_series_data(ac_ids),
  metadata_cols = "ac_amount"
)
out_ac <- summarize_alpha_gamma_ratios(
  div_ac,
  series = "ac",
  scenario_col = "ac_amount"
)


# Write output -----------------------------------------------------------

results <- bind_rows(out_core, out_hab, out_disp, out_ac)

out_file <- here("output/exploration/alpha_gamma_ratio_final.csv")
fwrite(results, out_file)
cat("\nAlpha/gamma ratio summary written to", out_file, "\n")
