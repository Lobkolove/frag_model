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
  end = vector("list", length(data_wide))
)
for (i in seq_along(data_wide)) {
  
  start <- Sys.time()

  sim_id <- data_wide[[i]]$sim_id[1]
  
  post_frag <- data_wide[[i]] %>%
    dplyr::filter(step_label == "post_fragmentation")
  
  end <- data_wide[[i]] %>%
    dplyr::filter(step_label == "final")
  
  fragmentation <- end$fragmentation[1]
  
  cat("Processing dataset", i, "( fragmentation =", fragmentation, "| ac_amount = " , data_wide[[i]]$ac_amount[1], ")\n")
  dd_results$post_frag[[i]] <- dist_decay(
    post_frag,
    dist_type = "euclidean",
    distvec = seq(0, 25, length.out = 200)
  )$smooth %>%
    dplyr::select(distance, similarity) %>%
    dplyr::mutate(sim_id = sim_id, fragmentation = fragmentation, step_label = "post_fragmentation")

  dd_results$end[[i]] <- dist_decay(
    end,
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



