library(tidyverse)
library(patchwork)
library(colorspace)

# Read in aggregated results
aggr_results <- readRDS("aggregated_state_results.rds")

# Averages and confidence intervals for plotting
sSBR_avg <- aggr_results %>%
  unnest(sSBR_main) %>%
  group_by(frag_level, samp_method, effort) %>%
  summarise(
    S_mean = mean(S), 
    S_low = quantile(S, 0.025),
    S_high = quantile(S, 0.975),
    .groups = "drop"
  ) %>%
  mutate(frag_level = factor(frag_level, levels = c("low", "high")),
         samp_method = factor(samp_method, labels = c("All habitat cells", "Chessboard sampling", "30 random habitat cells")))

dd_avg <- aggr_results %>%
  unnest(dd_main) %>%
  group_by(frag_level, samp_method, distance) %>%
  summarise(
    similarity_mean = mean(similarity),
    similarity_low = quantile(similarity, 0.025),
    similarity_high = quantile(similarity, 0.975),
    .groups = "drop"
  ) %>%
  mutate(frag_level = factor(frag_level, levels = c("low", "high")),
         samp_method = factor(samp_method, labels = c("All habitat cells", "Chessboard sampling", "30 random habitat cells")))


# Plotting
pal1 <- c(colorspace::lighten("midnightblue", 0.2), colorspace::lighten("violetred4", 0.2))

p1 <- ggplot(sSBR_avg, 
             aes(x = effort, y = S_mean, color = frag_level, fill = frag_level)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), alpha = 0.4, color = NA) +
  scale_color_manual(values = pal1) +
  scale_fill_manual(values = colorspace::lighten(pal1, 0.2)) +
  labs(x = "Cumulative area of convex hull", 
       y = "Cumulative species richness",
       title = "Spatially Constrained Rarefaction") +
  facet_wrap(~ samp_method, scales = "free", ncol = 3) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom")

p2 <- ggplot(dd_avg, 
             aes(x = distance, y = similarity_mean, color = frag_level, fill = frag_level)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = similarity_low, ymax = similarity_high), 
              alpha = 0.4, col = NA) +
  scale_color_manual(values = pal1) +
  scale_fill_manual(values = colorspace::lighten(pal1, 0.2)) +
  labs(x = "Euclidean distance", 
       y = "Similarity (1 - Bray-Curtis dissimilarity)",
       title = "Distance Decay of Similarity") +
  facet_wrap(~ samp_method, scales = "free", ncol = 3) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom")

# Combine: 2 rows x 3 columns = 6 panels
# Single legend shared across all panels
combined_plot <- p1 / p2 + 
  plot_layout(guides = "collect") & 
  guides(color = guide_legend(title = "Fragmentation", nrow = 2),
         fill = guide_legend(title = "Fragmentation", nrow = 2)) &
  theme(legend.position = "right")


print(combined_plot)

ggsave("aggr_state_combined_quantiles.png", 
       combined_plot, width = 10, height = 7.5, dpi = 300)
