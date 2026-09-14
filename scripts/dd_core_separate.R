# Distance decay of the core series as separate plots
#
# Reproduces the distance decay plots of the core series from
# scripts/main_analyses.R (sections 1.1.1 and 1.2.1), but draws one plot per
# time step instead of faceting by step_label.
# Everything up to the plotting is identical to main_analyses.R, so the curves
# are the same; only the layout differs.

library(here)
library(data.table)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tidyr)
library(vegan)

source("R/sim_select.R")
source("R/dist_decay.R")

# Plot settings
theme_set(theme_bw(base_size = 14))
pal_frag <- c("#4688ad", "#e9b14a", "#ac5384")

# Distances at which the fitted curves are evaluated
distvec <- seq(0, 25, length.out = 200)

# Read in log file
log_file <- here("output/simulations_log.csv")
log <- fread(log_file)


# 0. Helpers -------------------------------------------------------------

#' Read the sampled data of a set of simulations
#'
#' @param paths Character vector of file paths (relative to the project root).
#' @return A list of data tables, filtered to the two focal time steps.
read_core_data <- function(paths) {
  map(
    here(paths),
    ~ fread(.x) %>%
      dplyr::filter(step_label %in% c("post_fragmentation", "final"))
  )
}

#' Compute and summarise distance decay across replicates
#'
#' Computes the distance decay for every dataset and summarises the fitted
#' curves across replicates as mean and 95% quantile interval.
#'
#' @param data_list A list of data tables, as returned by `read_core_data()`.
#' @return A data frame with one row per fragmentation level, step and distance.
summarise_ddecay <- function(data_list) {
  map(
    data_list,
    ~ grouped_ddecay(
      model_sample = .x,
      distvec = distvec
    )
  ) |>
    bind_rows() |>
    dplyr::group_by(fragmentation, step_label, distance) |>
    dplyr::summarise(
      simi_low = quantile(similarity, 0.025, na.rm = TRUE),
      simi_high = quantile(similarity, 0.975, na.rm = TRUE),
      similarity = mean(similarity, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      fragmentation = factor(
        fragmentation,
        levels = c(0.2, 0.5, 0.8),
        labels = c("Low", "Medium", "High")
      ),
      step_label = factor(
        step_label,
        levels = c("post_fragmentation", "final"),
        labels = c("Post-fragmentation", "End of simulation")
      )
    )
}

#' Plot the distance decay of a single time step
#'
#' Draws the same curves as the faceted plots in clean_diversity.R, but for one
#' time step only. The y axis range is taken from the full data set so that the
#' separate plots stay comparable, exactly as the shared axis of a facet would.
#'
#' @param dd_data Summarised distance decay data (from `summarise_ddecay()`).
#' @param step Time step to plot; one of the levels of `step_label`.
#' @param title Plot title. Defaults to the time step.
#' @param ylim Common y axis range. Defaults to the range of `dd_data`.
#' @return A ggplot object.
plot_ddecay_step <- function(
  dd_data,
  step,
  title = step,
  ylim = range(c(dd_data$simi_low, dd_data$simi_high), na.rm = TRUE)
) {

  step <- match.arg(step, levels(dd_data$step_label))

  dd_step <- dd_data %>%
    dplyr::filter(step_label == step)

  ggplot(
    dd_step,
    aes(distance, similarity, color = fragmentation, fill = fragmentation)
  ) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(
      aes(ymin = simi_low, ymax = simi_high, fill = fragmentation),
      alpha = 0.2,
      color = NA
    ) +
    coord_cartesian(ylim = ylim) +
    labs(
      title = title,
      x = "Euclidean Distance",
      y = "Similarity (1 - Bray-Curtis dissimilarity)",
      color = "Level of\nFragmentation",
      fill = "Level of\nFragmentation"
    ) +
    scale_color_manual(values = pal_frag) +
    scale_fill_manual(values = pal_frag) +
    theme(legend.position = "right")
}


# 1. Core series ---------------------------------------------------------

# Assess core simulation IDs
core_ids <- sim_ids(habitat == 0.15, ac_amount == 0.7, dispersal_dist == 2)

# Check that no simulations share the same master seed
core_seeds <- log %>%
  filter(sim_id %in% core_ids) %>%
  pull(master_seed)
sum(duplicated(core_seeds))
# No duplicated master seeds, ready for analysis!


  ## 1.1 Full sample -------------------------------------------------------

# Select simulations, read in data and compute distance decay
paths_core_full <- sim_select(sim_id %in% core_ids, sampled = "all")
data_core_full <- read_core_data(paths_core_full)
dd_core_merged <- summarise_ddecay(data_core_full)

# Common y axis range, so the two separate plots remain comparable
ylim_full <- range(
  c(dd_core_merged$simi_low, dd_core_merged$simi_high),
  na.rm = TRUE
)

# One plot per time step
gg_dd_core_full_post <- plot_ddecay_step(
  dd_core_merged,
  step = "Post-fragmentation",
  ylim = ylim_full
)
gg_dd_core_full_post

gg_dd_core_full_final <- plot_ddecay_step(
  dd_core_merged,
  step = "End of simulation",
  ylim = ylim_full
)
gg_dd_core_full_final

ggsave(
  gg_dd_core_full_post,
  file = here("pics/dd_core_full_post_fragmentation.png"),
  width = 6,
  height = 4,
  dpi = 300
)

ggsave(
  gg_dd_core_full_final,
  file = here("pics/dd_core_full_final.png"),
  width = 6,
  height = 4,
  dpi = 300
)


  ## 1.2 Random samples ----------------------------------------------------

# Select simulations, read in data and compute distance decay
paths_core_random <- sim_select(sim_id %in% core_ids, sampled = "random")
data_core_random <- read_core_data(paths_core_random)
dd_core_random_merged <- summarise_ddecay(data_core_random)

# Common y axis range for the random samples
ylim_random <- range(
  c(dd_core_random_merged$simi_low, dd_core_random_merged$simi_high),
  na.rm = TRUE
)

# One plot per time step
gg_dd_core_random_post <- plot_ddecay_step(
  dd_core_random_merged,
  step = "Post-fragmentation",
  ylim = ylim_random
)
gg_dd_core_random_post

gg_dd_core_random_final <- plot_ddecay_step(
  dd_core_random_merged,
  step = "End of simulation",
  ylim = ylim_random
)
gg_dd_core_random_final

ggsave(
  gg_dd_core_random_post,
  file = here("pics/dd_core_random_post_fragmentation.png"),
  width = 6,
  height = 5.5,
  dpi = 300
)

ggsave(
  gg_dd_core_random_final,
  file = here("pics/dd_core_random_final.png"),
  width = 6,
  height = 5.5,
  dpi = 300
)
