# Fragmentation metrics of the core landscapes
#
# Summarises the geometry of the fragmented landscapes used in main_analyses.R
# (habitat = 0.15, ac = 0.7, dispersal_dist = 2; fragmentation = 0.2 / 0.5 / 0.8).
# Metrics are computed on the post-fragmentation grid of every core simulation
# and then summarised (mean, sd, range) per fragmentation level, so the numbers
# describe exactly the landscapes the diversity analyses are based on.
#
# The grid is toroidal, so patches are delineated with toroidal_clump() and all
# distances are toroidal. Edge cells are also defined on the torus (a habitat
# cell with at least one matrix cell among its four orthogonal neighbours); this
# differs from the cells$edges stored in the model states, which come from
# raster::boundaries() and treat the grid border as an edge.
#
# Outputs
#   output/frag_metrics_per_landscape.csv   one row per simulation
#   output/frag_metrics_summary.csv         one row per fragmentation level and metric
#   pics/frag_metrics.png                   distribution of every metric across replicates

library(here)
library(raster) # loaded first so that dplyr::select() masks raster::select()
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

source("R/sim_select.R")
source("R/toroidal_clump.R")
source("R/toroidal_dist.R")

# Plot settings
theme_set(theme_bw(base_size = 14))
pal_frag <- c( "#4688ad",  "#e9b14a", "#ac5384")

# Neighbourhood rule for patch delineation (4 = orthogonal neighbours only, as
# used by sample_cells() and the other scripts in this repository)
directions <- 4


# 1. Metric functions ----------------------------------------------------

# Shift a matrix by dr rows and dc columns with wrap-around, so that
# wrap_shift(m, 1, 0)[i, j] is m[i + 1, j] on the torus
wrap_shift <- function(m, dr, dc) {
  n <- nrow(m)
  k <- ncol(m)
  m[((seq_len(n) - 1 + dr) %% n) + 1, ((seq_len(k) - 1 + dc) %% k) + 1]
}

# Compute all metrics for one post-fragmentation grid (RasterLayer with NA for
# matrix cells). Returns a one-row data.frame.
landscape_metrics <- function(grid, directions = 4) {

  hab <- !is.na(raster::as.matrix(grid))
  n_hab <- sum(hab)
  n_total <- length(hab)
  grid_size <- nrow(hab)

  ## Edges (toroidal, orthogonal neighbours) ----
  # Number of habitat neighbours of every cell
  n_hab_nb <- wrap_shift(hab, 1, 0) + wrap_shift(hab, -1, 0) + wrap_shift(hab, 0, 1) + wrap_shift(hab, 0, -1)
  # Sides of habitat cells facing matrix
  matrix_sides <- (4 - n_hab_nb)[hab]

  n_edge <- sum(matrix_sides > 0)
  n_interior <- n_hab - n_edge
  edge_length <- sum(matrix_sides)

  ## Patches ----
  # Binary raster, so that habitat cells with an environmental value of exactly
  # 0 are not treated as background by raster::clump()
  hab_r <- grid
  hab_r[!is.na(grid)] <- 1
  patch_id <- raster::as.matrix(toroidal_clump(hab_r, directions = directions))[hab]
  patch_sizes <- as.numeric(table(patch_id))
  n_patches <- length(patch_sizes)

  ## Inter-patch distances (toroidal, between cell centres) ----
  # Distance between every pair of habitat cells; from there the minimum
  # distance between every pair of patches (edge-to-edge in cell units, where
  # 1 means the two patches are one cell apart with nothing in between; under
  # the 4-neighbour rule diagonally touching cells form separate patches at
  # distance sqrt(2))
  if (n_patches > 1) {
    coords <- which(hab, arr.ind = TRUE)
    d_cells <- as.matrix(toroidal_dist(coords, grid_size = grid_size))
    patches <- sort(unique(patch_id))
    d_patches <- sapply(patches, function(i) {
      # Minimum distance from patch i to every habitat cell, then to every patch
      d_i <- apply(d_cells[patch_id == i, , drop = FALSE], 2, min)
      tapply(d_i, patch_id, min)
    })
    diag(d_patches) <- NA
    enn <- apply(d_patches, 1, min, na.rm = TRUE)
    interpatch <- d_patches[upper.tri(d_patches)]
  } else {
    enn <- NA_real_
    interpatch <- NA_real_
  }

  data.frame(
    habitat_cells = n_hab,
    habitat_prop = n_hab / n_total,
    n_patches = n_patches,
    mean_patch_size = mean(patch_sizes),
    median_patch_size = median(patch_sizes),
    max_patch_size = max(patch_sizes),
    largest_patch_index = max(patch_sizes) / n_hab,
    prop_single_cell_patches = mean(patch_sizes == 1),
    enn_mean = mean(enn),
    interpatch_mean = mean(interpatch),
    edge_cells = n_edge,
    interior_cells = n_interior,
    prop_edge = n_edge / n_hab,
    edge_interior_ratio = ifelse(n_interior > 0, n_edge / n_interior, NA_real_),
    edge_length = edge_length,
    perimeter_area_ratio = edge_length / n_hab,
    like_adjacency = 1 - edge_length / (4 * n_hab)
  )
}

# Labels and definitions used in the summary table and the figure
metric_info <- tibble::tribble(
  ~metric,                     ~label,                              ~definition,
  "habitat_cells",             "Habitat cells",                     "Number of habitat cells",
  "habitat_prop",              "Habitat proportion",                "Habitat cells / all cells",
  "n_patches",                 "Number of patches",                 "Patches under the 4-neighbour rule on the torus",
  "mean_patch_size",           "Mean patch size (cells)",           "Mean number of cells per patch",
  "median_patch_size",         "Median patch size (cells)",         "Median number of cells per patch",
  "max_patch_size",            "Largest patch (cells)",             "Number of cells in the largest patch",
  "largest_patch_index",       "Largest patch index",               "Largest patch / all habitat cells",
  "prop_single_cell_patches",  "Single-cell patches (prop.)",       "Proportion of patches consisting of one cell",
  "enn_mean",                  "Mean nearest-neighbour distance",   "Mean over patches of the toroidal distance to the nearest other patch (cell centres, edge to edge)",
  "interpatch_mean",           "Mean inter-patch distance",         "Mean over all patch pairs of the minimum toroidal distance between them",
  "edge_cells",                "Edge cells",                        "Habitat cells with >= 1 orthogonal matrix neighbour (toroidal)",
  "interior_cells",            "Interior cells",                    "Habitat cells with 4 habitat neighbours (toroidal)",
  "prop_edge",                 "Edge cells (prop.)",                "Edge cells / habitat cells",
  "edge_interior_ratio",       "Edge-to-interior ratio",            "Edge cells / interior cells",
  "edge_length",               "Total edge length",                 "Number of habitat-matrix cell borders (toroidal)",
  "perimeter_area_ratio",      "Perimeter-area ratio",              "Total edge length / habitat cells",
  "like_adjacency",            "Like adjacency",                    "Proportion of habitat cell sides bordering habitat (1 = fully compact)"
)


# 2. Core landscapes -------------------------------------------------------

log <- fread(here("output/simulations_log.csv"))

# Same selection as in main_analyses.R
core_ids <- sim_ids(habitat == 0.15, ac_amount == 0.7, dispersal_dist == 2)
paths_core_state <- sim_select(sim_id %in% core_ids, file_type = "state")

# Read only the post-fragmentation grid of every state file
grids_core <- map(paths_core_state, function(path) {
  state <- readRDS(here(path))
  list(
    sim_id = state$post_fragmentation$meta$sim_id,
    fragmentation = state$post_fragmentation$meta$fragmentation,
    grid = state$post_fragmentation$grid
  )
})


# 3. Metrics per landscape ------------------------------------------------

metrics_landscape <- map_dfr(grids_core, function(ls) {
  cbind(
    data.frame(sim_id = ls$sim_id, fragmentation = ls$fragmentation),
    landscape_metrics(ls$grid, directions = directions)
  )
}) %>%
  arrange(fragmentation, sim_id)

fwrite(metrics_landscape, here("output/exploration/frag_metrics_per_landscape.csv"))


# 4. Summary per fragmentation level --------------------------------------

metrics_long <- metrics_landscape %>%
  pivot_longer(-c(sim_id, fragmentation), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metric_info$metric))

metrics_summary <- metrics_long %>%
  group_by(fragmentation, metric) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(metric_info, by = "metric") %>%
  select(fragmentation, metric, label, n, mean, sd, min, max, definition) %>%
  arrange(metric, fragmentation)

fwrite(metrics_summary, here("output/exploration/frag_metrics_summary.csv"))

# Wide overview: one row per metric, one column per fragmentation level (mean ± sd)
metrics_table <- metrics_summary %>%
  mutate(
    value = sprintf("%.2f ± %.2f", mean, sd),
    fragmentation = paste0("frag_", fragmentation)
  ) %>%
  select(metric, label, fragmentation, value) %>%
  pivot_wider(names_from = fragmentation, values_from = value)

print(as.data.frame(metrics_table), right = FALSE)


# 5. Figure ----------------------------------------------------------------

# Habitat amount is held constant and only listed as a check, so it is left out
# of the figure
plot_metrics <- c(
  "n_patches", "mean_patch_size", "median_patch_size", "largest_patch_index",
  "prop_single_cell_patches", "enn_mean", "interpatch_mean",
  "prop_edge", "edge_interior_ratio", "perimeter_area_ratio", "like_adjacency"
)

# NAs are landscapes with a single patch (no inter-patch distance) or without
# interior cells (no edge-to-interior ratio)
metrics_plot <- metrics_long %>%
  filter(metric %in% plot_metrics, !is.na(value)) %>%
  left_join(metric_info, by = "metric") %>%
  mutate(
    label = factor(label, levels = metric_info$label[metric_info$metric %in% plot_metrics]),
    fragmentation = factor(fragmentation, levels = c(0.2, 0.5, 0.8), labels = c("Low", "Medium", "High"))
  )

gg_frag_metrics <- ggplot(metrics_plot, aes(fragmentation, value, color = fragmentation)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_wrap(~ label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = pal_frag, guide = "none") +
  labs(x = "Level of fragmentation", y = "Metric value") +
  theme(strip.text = element_text(size = 10))
