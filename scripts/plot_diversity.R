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
theme_set(theme_bw(base_size = 14))

# Steps to analyse
steps <- c("post_fragmentation", "final")

# Fragmentation levels
frag_levels <- c("low" = 0.2, "medium" = 0.5, "high" = 0.8)


# Core series ------------------------------------------------------------

# Select simulation and get path for full sample file
paths_core_full <- sim_select(sampled = "all")

# Read in data
data_core_full <- map(here(paths_core_full), fread)

# Distance-decay results
# dd_results <- vector(mode = "list", length = length(steps))

div_core_test <- map(data_core_full, ~ compute_diversity(data = .x))

div_core_full <- vector(mode = "list", length = length(data_core_full))
# For each dataset:
for (i in seq_along(data_core_full)) {

  # Get current data and filter out rows with NA species_id
  data <- data_core_full[[i]] |> 
    filter(!is.na(species_id))

  # Diversity results
  div_results <- vector(mode = "list", length = length(steps))

  # For each step, filter and calculate alpha and gamma diversity
  for (j in seq_along(steps)) {

    # Filter for current step
    subset <- data |> 
      dplyr::filter(step_label == steps[j])

    # Assess species richness
    gamma <- tibble(
      scale = factor("landscape", levels = c("sample", "landscape")),
      richness = length(unique(subset[["species_id"]]))
    )
    alpha <- subset |> 
      group_by(sample_id) |> 
      summarise(richness = n_distinct(species_id)) |> 
      summarise(richness = mean(richness, na.rm = TRUE)) |> 
      mutate(scale = factor("sample", levels = c("sample", "landscape"))) |> 
      select(scale, richness)
    
    # Reformat to wide
    wide <- subset |> 
      tidyr::pivot_wider(
        names_from = species_id,
        values_from = n,
        values_fill = 0,
        names_prefix = "sp_",
        names_sort = TRUE
      )
    
    # Extract site x species table (remove sp_NA column only if present)
    spec_table <- wide |> 
      select(starts_with("sp_"))

    # Hill numbers
      ## alpha
    alpha <- alpha |> 
      mutate(
        hill_shannon = mean(as.numeric(renyi(spec_table, scales = 2, hill = TRUE))),
        hill_simpson = mean(as.numeric(renyi(spec_table, scales = 3, hill = TRUE))),
        evenness = log(hill_shannon) / log(richness)
      )
      ## gamma
    gamma <- gamma |> 
      mutate(
        hill_shannon = as.numeric(renyi(colSums(spec_table), scales = 2, hill = TRUE)),
        hill_simpson = as.numeric(renyi(colSums(spec_table), scales = 3, hill = TRUE)),
        evenness = hill_shannon / richness
      )
    
    # Store diversity results
    div_results[[j]] <- bind_rows(alpha, gamma) |>
      mutate(
        step_label = steps[j],
        sim_id = subset[["sim_id"]][1],
        fragmentation = factor(subset[["fragmentation"]][1], levels = frag_levels, labels = names(frag_levels)),
        step_label = factor(step_label, levels = steps, labels = c("Post-fragmentation", "End of simulation"))
      ) |> 
        select(sim_id, fragmentation, step_label, scale, richness, hill_shannon, hill_simpson, evenness) |> 
        pivot_longer(cols = richness:evenness, names_to = "index", values_to = "value")
                                            
  }

  # Merge and store results for this simulation (all steps)
  div_core_full[[i]] <- bind_rows(div_results)

}

diversity <- bind_rows(div_core_full)

# Plots

# Option 1: Boxplots for each diversity index, faceted by step and scale
div_plots <- vector(mode = "list", length = length(unique(diversity$index)))
names(div_plots) <- c("Species richness", "Hill Shannon", "Hill Simpson", "Evenness")
pal <- c("#332288", "#89a0c4", "#CC6677")

for (i in seq_along(div_plots)) {

  subset <- filter(diversity, index == unique(diversity$index)[i])

  div_plots[[i]] <- ggplot(subset, aes(x = fragmentation, y = value, color = fragmentation)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    facet_grid2(scale ~ step_label, scales = "free_y", independent = "y") +
    labs(x = "Level of Fragmentation", y = names(div_plots)[i]) +
    theme_bw(base_size = 14) +
    scale_color_manual(values = pal) +
    theme(legend.position = "none")

}

div_combined_plot <- wrap_plots(div_plots, ncol = 1, axes = "collect")
div_combined_plot

ggsave(div_combined_plot, filename = here("pics", "core_diversity.png"), width = 9, height = 9, dpi = 300)


# Option 2: Boxplots for each index, faceted by step and grouped by scale (local vs landscape)

div_plots2 <- vector(mode = "list", length = length(unique(diversity$scale)))
names(div_plots2) <- c("Local scale (alpha)", "Landscape scale (gamma)")

for (i in seq_along(div_plots2)) {

  subset <- filter(diversity, scale == unique(diversity$scale)[i]) |> 
    mutate(index = factor(index, levels = c("richness", "hill_shannon", "hill_simpson", "evenness"), 
                          labels = c("Species richness", "Hill Shannon", "Hill Simpson", "Evenness")))
  
  indices <- levels(subset$index)
  
  # subplots <- vector(mode = "list", length = length(unique(subset$index)))
  # for (j in seq_along(unique(subset$index))) {
    
  #   subplots[[j]] <- ggplot(filter(subset, index == indices[j]), aes(x = fragmentation, y = value, color = fragmentation)) +
  #     geom_boxplot() +
  #     geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  #     facet_wrap(~ step_label, scales = "free") +
  #     labs(x = "Level of Fragmentation", y = indices[j]) +
  #     theme_bw(base_size = 14) +
  #     scale_color_manual(values = pal) +
  #     theme(legend.position = "none")

  # }
  
  # div_plots2[[i]] <- wrap_plots(subplots, ncol = 1, axes = "collect", guides = "collect")

  div_plots2[[i]] <- ggplot(subset, aes(x = fragmentation, y = value, color = fragmentation)) +
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    facet_grid2(index ~ step_label, scales = "free_y", independent = "y", switch = "y") +
    labs(x = "Level of Fragmentation", y = names(div_plots2)[i]) +
    theme_bw(base_size = 14) +
    scale_color_manual(values = pal) +
    scale_y_continuous(
      name = "",
      sec.axis = dup_axis(
        name = names(div_plots2)[i],
        breaks = NULL,
        labels = NULL
      )
    ) +
    theme(
      # axis.title.y.right = element_text(angle = 90),
      legend.position = "none"
    )

}

div_plots2_combined <- wrap_plots(div_plots2, ncol = 1, axes = "collect_x")
div_plots2_combined


# Option 3: timestep as color, scale as facets

div_summary <- diversity |> 
  group_by(step_label, scale, fragmentation, index) |> 
  summarise(
    mean = mean(value, na.rm = TRUE),
    q_low = quantile(value, probs = 0.025, na.rm = TRUE),
    q_high = quantile(value, probs = 0.975, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  mutate(
    index = factor(index, levels = c("richness", "hill_shannon", "hill_simpson", "evenness"), 
                        labels = c("Species richness", "Hill Shannon", "Hill Simpson", "Evenness")),
    scale = factor(scale, levels = c("sample", "landscape"), labels = c("Sample scale", "Landscape scale"))
  )

gg_div_full_merged <- ggplot(div_summary, aes(fragmentation, mean, color = step_label)) +
  geom_pointrange(aes(ymin = q_low, ymax = q_high), alpha = 0.85) +
  facet_grid2(scale ~ index, scales = "free_y", independent = "y") +
  labs(title = "All habitat cells", x = "Level of Fragmentation", y = "Index value", color = "Time step") +
  scale_color_manual(values = pal[c(1,3)]) +
  theme(legend.position = "bottom")
gg_div_full_merged

ggsave(gg_div_full_merged, filename = here("pics", "diversity_merged_fullsamp.png"), width = 10, height = 8, dpi = 300)

# Random sampling



# dispersal distance -----------------------------------------------------

# Select simulations with varying dispersal distances
paths <- sim_select(dispersal_dist != 2, sampled = "all")

# Read in data
data_list <- map(here(paths), fread)

# Dispersal distances
dispersal_dists <- c("short" = 1, "medium" = 4, "long" = 8)

# Initialize list to store diversity results
div_list <- vector(mode = "list", length = length(data_list))

for (i in seq_along(data_list)) {

  # Get current data and filter out rows with NA species_id
  data <- data_list[[i]] |> 
    filter(!is.na(species_id))

  # Diversity results
  div_results <- vector(mode = "list", length = length(steps))

  # For each step, filter and calculate alpha and gamma diversity
  for (j in seq_along(steps)) {

    # Filter for current step
    subset <- data |> 
      dplyr::filter(step_label == steps[j])

    # Assess species richness
    gamma <- tibble(
      scale = factor("landscape", levels = c("sample", "landscape")),
      richness = length(unique(subset[["species_id"]]))
    )
    alpha <- subset |> 
      group_by(sample_id) |> 
      summarise(richness = n_distinct(species_id)) |> 
      summarise(richness = mean(richness, na.rm = TRUE)) |> 
      mutate(scale = factor("sample", levels = c("sample", "landscape"))) |> 
      select(scale, richness)
    
    div_results[[j]] <- bind_rows(alpha, gamma) |>
      mutate(
        step_label = steps[j],
        sim_id = subset[["sim_id"]][1],
        fragmentation = factor(subset[["fragmentation"]][1], levels = frag_levels, labels = names(frag_levels)),
        dispersal_dist = factor(subset[["disp_dist"]][1], levels = dispersal_dists, labels = names(dispersal_dists)),
        step_label = factor(step_label, levels = steps, labels = c("Post-fragmentation", "End of simulation"))
      ) |> 
        select(sim_id, step_label, scale, fragmentation, dispersal_dist, richness)

  }

  div_list[[i]] <- bind_rows(div_results)

}

diversity <- bind_rows(div_list) %>% 
  group_by(step_label, scale, fragmentation, dispersal_dist) %>% 
  summarise(
    S_low = quantile(richness, probs = 0.025),
    S_high = quantile(richness, probs = 0.975),
    richness = mean(richness),
  ) %>% 
  ungroup()

ggplot(diversity, aes(x = fragmentation, y = richness, color = step_label)) +
  geom_pointrange(aes(ymin = S_low, ymax = S_high)) +
  # geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  facet_grid(scale ~ dispersal_dist, scales = "free_y") +
  labs(x = "Level of Fragmentation", y = "Species richness") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none")




