# Repeat the analysis of the "1.2 Random samples" subsection of
# clean_diversity.R, but draw the 30 cells from the full sample data so that the
# *same* cells are followed at both focal time steps.
#
# The pre-drawn "_rand" files cannot do this: sample_cells() calls
# sample(habitat_cells, n_samples) separately for every recorded state, so the
# post-fragmentation and final samples of a simulation are independent draws
# (they typically share only a handful of cells). Sampling from the "_all" files
# instead keeps cell identity fixed, which turns the comparison between the two
# time steps into a paired one.

library(here)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(patchwork)
library(stringr)
library(tidyr)
library(purrr)
library(vegan)

source("R/sim_select.R")
source("R/compute_diversity.R")
source("R/dist_decay.R")

# Plot settings
theme_set(theme_bw(base_size = 14))
pal_frag <- c("#4688ad", "#e9b14a", "#ac5384")
pal_geodem <- c("#503f80", "#dd7e8d")

# Sampling settings
n_samples <- 30
sample_seed <- 1312

# Focal time steps
steps <- c("post_fragmentation", "final")

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

# Read in data
data_core_full <- map(
  here(paths_core_full),
  ~ fread(.x) %>%
    dplyr::filter(step_label %in% steps)
)


  ## 1.1 Fixed random subsample --------------------------------------------

#' Draw a random subset of cells that is shared across time steps
#'
#' Restricts a simulation to `n_samples` habitat cells drawn at random from the
#' cells that are present at every focal step, so that the same cells are
#' compared across steps. Cells are kept regardless of whether they are
#' occupied, matching the behaviour of sample_cells(method = "random").
#'
#' @param data A sampled-data table of a single simulation ("_all" file).
#' @param n_samples Number of cells to draw.
#' @param steps Focal time steps the drawn cells have to be present at.
#' @param seed Seed for the draw, for reproducibility.
#'
#' @return `data`, filtered to the drawn cells.
subsample_cells <- function(data, n_samples = 30, steps = c("post_fragmentation", "final"), seed = NULL) {

  data <- data %>%
    dplyr::filter(step_label %in% steps)

  # Cells recorded at every focal step
  shared_cells <- data %>%
    dplyr::distinct(step_label, cell_id) %>%
    dplyr::count(cell_id) %>%
    dplyr::filter(n == length(steps)) %>%
    dplyr::pull(cell_id)

  if (length(shared_cells) < n_samples) {
    stop(
      "Only ", length(shared_cells), " cells are shared across all focal steps, ",
      "but ", n_samples, " were requested."
    )
  }

  if (!is.null(seed)) set.seed(seed)
  drawn <- sample(shared_cells, n_samples)

  data %>%
    dplyr::filter(cell_id %in% drawn)
}

# Draw the same 30 cells per simulation for both focal time steps.
# The seed is offset by sim_id so every simulation gets its own reproducible draw.
data_core_fixed <- map(
  data_core_full,
  ~ subsample_cells(
    .x,
    n_samples = n_samples,
    steps = steps,
    seed = sample_seed + .x$sim_id[1]
  )
)

# Sanity check: the same cells at both steps, in every simulation
stopifnot(
  all(map_lgl(
    data_core_fixed,
    ~ {
      cells_by_step <- split(.x$cell_id, .x$step_label)
      length(unique(cells_by_step[[1]])) == n_samples &&
        setequal(cells_by_step[[1]], cells_by_step[[2]])
    }
  ))
)


  ## 1.2 Diversity indices -------------------------------------------------

# Compute diversity indices for each dataset in the core series.
# cell_id is used as the sample identifier so that alpha diversity is
# aggregated over the drawn cells.
div_core_fixed <- map(
  data_core_fixed,
  ~ compute_diversity(data = .x, sample_col = "cell_id")
)

# Merge results into a single data frame
div_core_fixed <- bind_rows(div_core_fixed)

# Create summary data frame for plotting
div_core_fixed_summary <- div_core_fixed |>
  group_by(step_label, scale, fragmentation, index) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    q_low = quantile(value, probs = 0.025, na.rm = TRUE),
    q_high = quantile(value, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Pointrange plot with time step by color and facet grid by scale and index
gg_div_core_fixed <- ggplot(div_core_fixed_summary, aes(fragmentation, mean, color = step_label)) +
  geom_pointrange(aes(ymin = q_low, ymax = q_high), alpha = 0.85) +
  facet_grid2(scale ~ index, scales = "free_y", independent = "y") +
  labs(x = "Level of Fragmentation", y = "Index value", color = "Time step") +
  scale_color_manual(values = pal_geodem) +
  theme(legend.position = "bottom")
gg_div_core_fixed


  ## 1.3 Distance decay ----------------------------------------------------

# Compute distance decay for each dataset in the core series
dd_core_fixed <- map(
  data_core_fixed,
  ~ grouped_ddecay(
    model_sample = .x,
    distvec = seq(0, 25, length.out = 200)
  )
)

# Merge results into a single data frame
dd_core_fixed_merged <- bind_rows(dd_core_fixed) |>
  dplyr::group_by(fragmentation, step_label, distance) |>
  dplyr::summarise(
    simi_low = quantile(similarity, 0.025, na.rm = TRUE),
    simi_high = quantile(similarity, 0.975, na.rm = TRUE),
    similarity = mean(similarity, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    fragmentation = factor(fragmentation, levels = c(0.2, 0.5, 0.8), labels = c("Low", "Medium", "High")),
    step_label = factor(step_label, levels = c("post_fragmentation", "final"), labels = c("Post-fragmentation", "End of simulation"))
  )

# Plot distance decay curves
gg_dd_core_fixed <- ggplot(
  dd_core_fixed_merged,
  aes(distance, similarity, color = fragmentation, fill = fragmentation)
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = simi_low, ymax = simi_high, fill = fragmentation),
    alpha = 0.2,
    color = NA
  ) +
  facet_wrap(~step_label) +
  labs(
    x = "Euclidean Distance",
    y = "Similarity (1 - Bray-Curtis dissimilarity)",
    color = "Level of Fragmentation",
    fill = "Level of Fragmentation"
  ) +
  scale_color_manual(values = pal_frag) +
  scale_fill_manual(values = pal_frag) +
  theme(legend.position = "bottom")
gg_dd_core_fixed


  ## 1.4 Combined plot -----------------------------------------------------

# Diversity indices on top (A), distance decay below (B).
# Note that plot_annotation() has to be added with `+` rather than `&`: combined
# with free() the `&` form errors out in patchwork ("first argument must be a
# vector"), since `&` distributes the annotation over the subplots.
gg_core_fixed_combined <- (free(gg_div_core_fixed) / gg_dd_core_fixed) +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(title = "30 random habitat cells, fixed across time steps", tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))
gg_core_fixed_combined

ggsave(
  gg_core_fixed_combined,
  file = here("pics/core_fixed_cells_combined.png"),
  width = 10,
  height = 10,
  dpi = 300
)
