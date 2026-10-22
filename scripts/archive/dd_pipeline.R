library(tidyverse)
source("R/dist_decay.R")
source("R/toroidal_dist.R")

# Read in all full sampled datasets with an autocorrelation of 0.7
sampled_files <- list.files("output/sampled_data/", 
                            pattern = "ac0\\.7.*samp_cb\\.csv$", 
                            full.names = TRUE)

# Create a list with all datasets
data_list <- map(sampled_files, read_csv)

# Convert to wide format for spatial curve analyses
data_wide <- data_list %>%
  map(
    ~ .x %>%
      tidyr::pivot_wider(
        names_from = species_id,
        values_from = n,
        values_fill = 0,
        names_prefix = "sp_",
        names_sort = TRUE
      )
  )

# Run distance decay for each dataset
dd_results <- list(
  post_frag = vector("list", length(data_wide)),
  final = vector("list", length(data_wide))
)
for (i in seq_along(data_wide)) {
  
  start <- Sys.time()

  sim_id <- data_wide[[i]]$sim_id[1]
  
  post_frag <- data_wide[[i]] %>%
    dplyr::filter(step_label == "post_fragmentation")
  
  final <- data_wide[[i]] %>%
    dplyr::filter(step_label == "final")
  
  fragmentation <- final$fragmentation[1]
  
  cat("Processing dataset", i, "( fragmentation =", fragmentation, "| ac_amount = " , data_wide[[i]]$ac_amount[1], ")\n")
  dd_results$post_frag[[i]] <- dist_decay(
    post_frag,
    dist_type = "euclidean",
    distvec = seq(0, 25, length.out = 200)
  )$smooth %>%
    dplyr::select(distance, similarity) %>%
    dplyr::mutate(sim_id = sim_id, fragmentation = fragmentation, step_label = "post_fragmentation")

  dd_results$final[[i]] <- dist_decay(
    final,
    dist_type = "euclidean",
    distvec = seq(0, 25, length.out = 200)
  )$smooth %>%
    dplyr::select(distance, similarity) %>%
    dplyr::mutate(sim_id = sim_id, fragmentation = fragmentation, step_label = "final")
  
  end <- Sys.time()
  cat("Time taken:", round(difftime(end, start, units = "secs"), 2), "seconds\n")
}

# Average distance decay curves across time step for each fragmentation level
dd_summaries <- map(dd_results, ~ .x %>%
  bind_rows() %>%
  group_by(step_label, fragmentation, distance) %>%
  summarise(simi_low = quantile(similarity, 0.025, na.rm = TRUE),
            simi_high = quantile(similarity, 0.975, na.rm = TRUE),
            similarity = mean(similarity, na.rm = TRUE), .groups = "drop") %>%
  mutate(fragmentation = factor(fragmentation, levels = c(0.2, 0.5, 0.8), labels = c("Low", "Medium", "High"))))

dd_summary <- bind_rows(dd_summaries) %>% 
  mutate(step_label = factor(step_label, levels = c("post_fragmentation", "final"), labels = c("Post-fragmentation", "End of simulation")))

# Plot distance decay curves
pal <- c("midnightblue", "darkgoldenrod", "violetred4")
ggplot(dd_summary, aes(x = distance, y = similarity, color = fragmentation, fill = fragmentation)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_ribbon(aes(ymin = simi_low, ymax = simi_high), alpha = 0.2, color = NA) +
  facet_wrap(~ step_label) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = "Euclidean Distance", y = "Similarity") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")



# Assess alpha and gamma div for each dataset
SR_results <- list(
  post_frag = vector("list", length(data_list)),
  final = vector("list", length(data_list))
)
for (i in seq_along(data_list)) {

  post_frag <- data_list[[i]] %>% 
    dplyr::filter(step_label == "post_fragmentation")

  pf_summary <- post_frag %>% 
    dplyr::group_by(sim_id, fragmentation, sample_id) %>% 
    dplyr::summarise(richness = n(), .groups = "drop")
  
  pf_alpha <- mean(pf_summary$richness, na.rm = TRUE)
  pf_gamma <- length(unique(post_frag$species_id))
  
  SR_results$post_frag[[i]] <- list(summary = pf_summary,
                                    alpha = pf_alpha,
                                    gamma = pf_gamma)
  
  final <- data_list[[i]] %>% 
      dplyr::filter(step_label == "final")
  
  fin_summary <- final %>% 
    dplyr::group_by(sim_id, fragmentation, sample_id) %>% 
    dplyr::summarise(richness = n(), .groups = "drop")

  fin_alpha <- mean(fin_summary$richness, na.rm = TRUE)
  fin_gamma <- length(unique(final$species_id))

  sim_id <- first(final$sim_id)
  fragmentation <- first(final$fragmentation)

  SR_results$final[[i]] <- list(sim_id = sim_id,
                                fragmentation = fragmentation,
                                summary = fin_summary,
                                alpha = fin_alpha,
                                gamma = fin_gamma)

}

# Aggregate species richness data into tidy format for plotting
# Alpha diversity (local scale) - per sample
alpha_tidy <- bind_rows(
  # Post-fragmentation
  SR_results$post_frag %>%
    map_df(
      ~ .x$summary %>%
        mutate(step_label = "post_fragmentation"),
      .id = "dataset_idx"
    ),
  
  # Final
  SR_results$final %>%
    map_df(
      ~ .x$summary %>%
        mutate(step_label = "final"),
      .id = "dataset_idx"
    )
) %>%
  mutate(scale = "Local (alpha)",
         richness_value = richness) %>%
  select(sim_id, fragmentation, step_label, scale, richness_value)

# Gamma diversity (landscape scale) - per dataset
gamma_tidy <- bind_rows(
  # Post-fragmentation
  SR_results$post_frag %>%
    map_df(
      ~ {
        tibble(
          sim_id = .x$summary$sim_id[1],
          fragmentation = .x$summary$fragmentation[1],
          step_label = "post_fragmentation",
          scale = "Landscape (gamma)",
          richness_value = .x$gamma
        )
      }
    ),
  
  # Final
  SR_results$final %>%
    map_df(
      ~ {
        tibble(
          sim_id = .x$sim_id,
          fragmentation = .x$summary$fragmentation[1],
          step_label = "final",
          scale = "Landscape (gamma)",
          richness_value = .x$gamma
        )
      }
    )
)

# Combine alpha and gamma
sr_combined <- bind_rows(alpha_tidy, gamma_tidy) %>%
  mutate(fragmentation = factor(fragmentation, levels = c(0.2, 0.5, 0.8), labels = c("Low", "Medium", "High")),
         step_label = factor(step_label, levels = c("post_fragmentation", "final"), labels = c("Post-fragmentation", "End of simulation")),
         scale = factor(scale, levels = c("Local (alpha)", "Landscape (gamma)")))

# Plot species richness with violin plots
ggplot(sr_combined, aes(x = fragmentation, y = richness_value, color = fragmentation, fill = fragmentation)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  facet_grid(scale ~ step_label, scales = "free_y") +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = "Fragmentation Level", y = "Species Richness") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
