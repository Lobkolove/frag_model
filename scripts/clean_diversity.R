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
pal_frag <- c( "#4688ad",  "#e9b14a", "#ac5384")
pal_geodem <- c("#503f80", "#dd7e8d")

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


  ## 1.1 Full sample -------------------------------------------------------------

# Select simulation and get path for full sample file
paths_core_full <- sim_select(sim_id %in% core_ids, sampled = "all")

# Read in data
data_core_full <- map(
  here(paths_core_full), 
  ~ fread(.x) %>%
    dplyr::filter(step_label %in% c("post_fragmentation", "final"))
)

    ### 1.1.1 Distance decay -----

# Compute distance decay for each dataset in the core series
dd_core_full <- map(
  data_core_full, 
  ~ grouped_ddecay(
    model_sample = .x,
    distvec = seq(0, 25, length.out = 200)
  )
)

# Merge results into a single data frame
dd_core_merged <- bind_rows(dd_core_full) |>
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

      #### 1.1.1.1 Fragmentation as color, step as facet -----

gg_dd_core_full <- ggplot(
  dd_core_merged,
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
gg_dd_core_full

    ### 1.1.2 Diversity indices -----
  
# Compute diversity indices for each dataset in the core series
div_core_full <- map(data_core_full, ~ compute_diversity(data = .x))

# Merge results into a single data frame
div_core_full <- bind_rows(div_core_full)


      #### 1.1.2.1 Pointrange plot, steps merged -----

# Create summary data frame for plotting
div_core_full_summary <- div_core_full |> 
  group_by(step_label, scale, fragmentation, index) |> 
  summarise(
    mean = mean(value, na.rm = TRUE),
    q_low = quantile(value, probs = 0.025, na.rm = TRUE),
    q_high = quantile(value, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Create pointrange plot with geodem by color and facet grid by scale and index
gg_div_core_full_merged <- ggplot(div_core_full_summary, aes(fragmentation, mean, color = step_label)) +
  geom_pointrange(aes(ymin = q_low, ymax = q_high), alpha = 0.85) +
  facet_grid2(scale ~ index, scales = "free_y", independent = "y") +
  labs(x = "Level of Fragmentation", y = "Index value", color = "Time step") +
  scale_color_manual(values = pal_geodem) +
  theme(legend.position = "bottom")
gg_div_core_full_merged

    ### 1.1.3 Combined plots -----

gg_core_combined <- (gg_dd_core_full / free(gg_div_core_full_merged)) +
  plot_layout(heights = c(1, 1.5)) &
  plot_annotation(title = "All habitat cells", tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))
gg_core_combined

ggsave(gg_core_combined, file = here("pics/core_combined.png"), width = 10, height = 10, dpi = 300)

  ## 1.2 Random samples -------------------------------------------------------------

# Select simulation and get path for random sample file
paths_core_random <- sim_select(sim_id %in% core_ids, sampled = "random")

# Read in data
data_core_random <- map(here(paths_core_random), fread)

# Filter to only include post-fragmentation and final steps
data_core_random <- map(data_core_random, ~ dplyr::filter(.x, step_label %in% c("post_fragmentation", "final")))

    ### 1.2.1 Distance decay -----

# Compute distance decay for each dataset in the core series
dd_core_random <- map(
  data_core_random, 
  ~ grouped_ddecay(
    model_sample = .x,
    distvec = seq(0, 25, length.out = 200)
  )
)

# Merge results into a single data frame
dd_core_random_merged <- bind_rows(dd_core_random) |>
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

# Plot distance decay curves for random samples
gg_dd_core_random <- ggplot(
  dd_core_random_merged,
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
gg_dd_core_random

    ### 1.2.2 Diversity indices -----

# Compute diversity indices for each dataset in the core series
div_core_random <- map(data_core_random, ~ compute_diversity(data = .x))

# Merge results into a single data frame
div_core_random <- bind_rows(div_core_random)

# Create summary data frame for plotting
div_core_random_summary <- div_core_random |>
  group_by(step_label, scale, fragmentation, index) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    q_low = quantile(value, probs = 0.025, na.rm = TRUE),
    q_high = quantile(value, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Pointrange plot for random samples, with geodem by color and facet grid by scale and index
gg_div_core_random_merged <- ggplot(div_core_random_summary, aes(fragmentation, mean, color = step_label)) +
  geom_pointrange(aes(ymin = q_low, ymax = q_high)) +
  facet_grid2(scale ~ index, scales = "free_y", independent = "y") +
  labs(x = "Level of Fragmentation", y = "Index value", color = "Time step") +
  scale_color_manual(values = pal_geodem) +
  theme(legend.position = "bottom")
gg_div_core_random_merged

# Combine the two plots into a single figure with patchwork, adding plot annotations (A and B)
# and only one x axis label for the combined figure
gg_div_core_random_combined <- (gg_div_core_full_merged / gg_div_core_random_merged) +
  plot_layout(axis_titles = "collect_x") &
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14), legend.position = "bottom")
gg_div_core_random_combined

ggsave(
  gg_div_core_random_combined,
  file = here("pics/diversity_core_random_combined.png"),
  width = 10,
  height = 8,
  dpi = 300
)

    ### 1.2.3 Combined plots for random samples -----

gg_core_random_combined <- (gg_dd_core_random / free(gg_div_core_random_merged)) +
  plot_layout(heights = c(1, 1.5)) &
  plot_annotation(title = "30 random habitat cells", tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))
gg_core_random_combined

ggsave(
  gg_core_random_combined,
  file = here("pics/core_random_combined.png"),
  width = 10,
  height = 10,
  dpi = 300
)


# 2. Habitat amount series ---------------------------------------------------------

# IDs
habitat_ids <- sim_ids(habitat %in% c(0.1, 0.3, 0.5))

# Check master seeds for duplicates
habitat_seeds <- log %>%
  dplyr::filter(sim_id %in% habitat_ids) %>%
  dplyr::pull(master_seed)
sum(duplicated(habitat_seeds))
# No seed duplicates, ready for analysis!

# Get paths for full samples
paths_hab <- sim_select(sim_id %in% habitat_ids, sampled = "all")

# Read in data
data_hab <- map(here(paths_hab), fread)

# Filter to only include post-fragmentation and final steps
data_hab <- map(data_hab, ~ dplyr::filter(.x, step_label %in% c("post_fragmentation", "final")))

  ## 2.1 Distance decay ---------------------------------------------------------

# Compute distance decay for each dataset in the habitat amount series
dd_hab <- map(
  data_hab, 
  ~ grouped_ddecay(
    model_sample = .x,
    group_cols = c("fragmentation", "habitat", "step_label"),
    distvec = seq(0, 25, length.out = 200)
  )
)

# Merge results into a single data frame
dd_hab_merged <- bind_rows(dd_hab) |>
  dplyr::group_by(fragmentation, habitat, step_label, distance) |>
  dplyr::summarise(
    simi_low = quantile(similarity, 0.025, na.rm = TRUE),
    simi_high = quantile(similarity, 0.975, na.rm = TRUE),
    similarity = mean(similarity, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    fragmentation = factor(fragmentation, levels = c(0.2, 0.5, 0.8), labels = c("Low", "Medium", "High")),
    habitat = factor(habitat, levels = c(0.1, 0.3, 0.5), labels = c("10% habitat", "30% habitat", "50% habitat")),
    step_label = factor(step_label, levels = c("post_fragmentation", "final"), labels = c("Post-fragmentation", "End of simulation"))
  )

# Plot distance decay curves for habitat amount series
gg_dd_hab <- ggplot(dd_hab_merged, aes(x = distance, y = similarity, color = fragmentation)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = simi_low, ymax = simi_high, fill = fragmentation), alpha = 0.2, color = NA) +
  facet_grid2(step_label ~ habitat, scales = "free_y") +
  labs(
    x = "Euclidean Distance",
    y = "Similarity (1 - Bray-Curtis dissimilarity)",
    color = "Level of\nfragmentation",
    fill = "Level of\nfragmentation"
  ) +
  scale_color_manual(values = pal_frag) +
  scale_fill_manual(values = pal_frag)
gg_dd_hab

  ## 2.2 Diversity indices ---------------------------------------------------------

# Compute diversity indices for each dataset in the habitat amount series
div_hab <- map(data_hab, ~ compute_diversity(data = .x, metadata_cols = "habitat"))

# Merge results into a single data frame
div_hab <- bind_rows(div_hab)

# Create summary data frame for plotting
div_hab_summary <- div_hab |>
  filter(index == "richness") |>
  group_by(step_label, scale, fragmentation, habitat) |>
  summarise(
    richness = mean(value, na.rm = TRUE),
    S_low = quantile(value, probs = 0.025, na.rm = TRUE),
    S_high = quantile(value, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  mutate(
    habitat = factor(habitat, levels = c(0.1, 0.3, 0.5), labels = c("10% habitat", "30% habitat", "50% habitat")),
  )

# Pointrange plot for habitat amount series, with geodem by color and facet grid by scale and habitat
gg_div_hab <- ggplot(div_hab_summary, aes(x = fragmentation, y = richness, color = step_label)) +
  geom_pointrange(aes(ymin = S_low, ymax = S_high), alpha = 0.85) +
  facet_grid2(scale ~ habitat, scales = "free_y") +
  labs(x = "Level of fragmentation", y = "Species richness", color = "Time step") +
  scale_color_manual(values = pal_geodem)
gg_div_hab

  ## 2.3 Combined plots for habitat amount series -----

gg_hab_combined <- free(gg_dd_hab, side = "l") / free(gg_div_hab, side = "l") +
  # plot_layout(guides = "keep") &
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14), legend.position = "right")
gg_hab_combined  

ggsave(
  gg_hab_combined,
  file = here("pics/habitat_combined_aligned.png"),
  width = 10,
  height = 10,
  dpi = 300
)


# 3. Dispersal distance series ---------------------------------------------------------

# IDs
dispersal_ids <- sim_ids(dispersal_dist %in% c(1, 4, 8))

# Check master seeds for duplicates
dispersal_seeds <- log %>%
  dplyr::filter(sim_id %in% dispersal_ids) %>%
  dplyr::pull(master_seed)
sum(duplicated(dispersal_seeds))
# No seed duplicates, ready for analysis!

# Get paths for full samples
dispersal_paths <- sim_select(sim_id %in% dispersal_ids, sampled = "all")

# Read in data
dispersal_data <- map(here(dispersal_paths), fread)

# Filter to only include post-fragmentation and final steps
dispersal_data <- map(dispersal_data, ~ dplyr::filter(.x, step_label %in% c("post_fragmentation", "final")))

  ## 3.1 Distance decay ---------------------------------------------------------

# Compute distance decay for each dataset in the dispersal distance series
dd_dispersal <- purrr::map(
  seq_along(dispersal_data),
  \(i) {

    tryCatch(
      grouped_ddecay(
        model_sample = dispersal_data[[i]],
        group_cols = c(
          "fragmentation",
          "disp_dist",
          "step_label"
        ),
        distvec = seq(0, 25, length.out = 200)
      ),
      error = function(e) {

        message(
          "Skipping dataset ", i, ": ",
          conditionMessage(e)
        )

        NULL
      }
    )
  }
)

skipped_datasets <- which(
  purrr::map_lgl(dd_dispersal, is.null)
)

n_skipped <- length(skipped_datasets)

# One dataset was skipped due to an error in the distance decay computation.
# The skipped dataset is 54, so we will exclude the whole replicate (sims 46-54) from the analysis.
dd_dispersal[46:54] <- NULL

# Merge results into a single data frame
dd_dispersal_merged <- bind_rows(dd_dispersal) |>
  dplyr::group_by(fragmentation, disp_dist, step_label, distance) |>
  dplyr::summarise(
    simi_low = quantile(similarity, 0.025, na.rm = TRUE),
    simi_high = quantile(similarity, 0.975, na.rm = TRUE),
    similarity = mean(similarity, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    fragmentation = factor(fragmentation, levels = c(0.2, 0.5, 0.8), labels = c("Low", "Medium", "High")),
    disp_dist = factor(disp_dist, levels = c(1, 4, 8), labels = c("1 cell", "4 cells", "8 cells")),
    step_label = factor(step_label, levels = c("post_fragmentation", "final"), labels = c("Post-fragmentation", "End of simulation"))
  )

# Plot distance decay curves for dispersal distance series
gg_dd_dispersal <- ggplot(dd_dispersal_merged, aes(x = distance, y = similarity, color = fragmentation)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = simi_low, ymax = simi_high, fill = fragmentation), alpha = 0.2, color = NA) +
  facet_grid2(step_label ~ disp_dist, scales = "free_y") +
  labs(
    x = "Euclidean Distance",
    y = "Similarity (1 - Bray-Curtis dissimilarity)",
    color = "Level of\nfragmentation",
    fill = "Level of\nfragmentation"
  ) +
  scale_color_manual(values = pal_frag) +
  scale_fill_manual(values = pal_frag)
gg_dd_dispersal

  ## 3.2 Diversity indices ---------------------------------------------------------

# Compute diversity indices for each dataset in the dispersal distance series
div_dispersal <- map(dispersal_data, ~ compute_diversity(data = .x, metadata_cols = "disp_dist"))

# Merge results into a single data frame
div_dispersal <- bind_rows(div_dispersal)

# Create summary data frame for plotting
div_dispersal_summary <- div_dispersal |>
  filter(index == "richness") |>
  group_by(step_label, scale, fragmentation, disp_dist) |>
  summarise(
    richness = mean(value, na.rm = TRUE),
    S_low = quantile(value, probs = 0.025, na.rm = TRUE),
    S_high = quantile(value, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  mutate(
    disp_dist = factor(disp_dist, levels = c(1, 4, 8), labels = c("1 cell", "4 cells", "8 cells")),
  )

# Pointrange plot for dispersal distance series, with geodem by color and facet grid by scale and dispersal distance
gg_div_dispersal <- ggplot(div_dispersal_summary, aes(x = fragmentation, y = richness, color = step_label)) +
  geom_pointrange(aes(ymin = S_low, ymax = S_high), alpha = 0.85) +
  facet_grid2(scale ~ disp_dist, scales = "free_y") +
  labs(x = "Level of fragmentation", y = "Species richness", color = "Time step") +
  scale_color_manual(values = pal_geodem)
gg_div_dispersal

  ## 3.3 Combined plots for dispersal distance series -----

gg_dispersal_combined <- free(gg_dd_dispersal, side = "l") / free(gg_div_dispersal, side = "l") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14), legend.position = "right")
gg_dispersal_combined
