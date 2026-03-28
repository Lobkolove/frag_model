library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(purrr)
library(colorspace)
library(patchwork)
library(raster)
library(readr)
source("Model/src/landscape.R")
source("Model/src/fragmentation.R")
source("R/sample_cells.R")
source("R/toroidal_clump.R")
source("R/dist_decay.R")
source("R/sSBR.R")

habitat <- 0.15
frag_levels <- c(low = 0.25, high = 0.75)
samp_methods <- c("all", "chessboard", "random")
n_reps <- 30

# Load states
aggr_state <- readRDS("data-raw/states_frag_0.5_sim_7.rds")$pre_fragmentation
rand_state <- readRDS("data-raw/random_state_large.rds")$start

# Function for single iteration
process_single_rep <- function(state, frag_level, samp_method, rep_id) {
  # Set unique seed for this replication
  seed_fragment <- 100 + rep_id
  
  # Fragment
  frag_state <- fragment(full_state = state, 
                         habitat = habitat, 
                         fragmentation = frag_levels[[frag_level]],
                         seed_fragment = seed_fragment)
  
  # Sample with explicit method
  sampled <- sample_cells(full_state = frag_state, method = samp_method, n_samples = 30)
  
  # Extract main curves only
  sSBR_main <- sSBR(model_sample = sampled, method = "area", 
                    distvec = seq(0, 625, length = 200))$smooth %>% 
    dplyr::select(effort, S)
  
  dd_main <- dist_decay(model_sample = sampled, dist_type = "euclidean", 
                        distvec = seq(0, 25, length = 200))$smooth %>% 
    dplyr::select(distance, similarity)
  
  # Return tidy row
  tibble(
    frag_level = frag_level,
    samp_method = samp_method,
    rep = rep_id,
    sSBR_main = list(sSBR_main),
    dd_main = list(dd_main)
  )
}

# Create parameter grid for iterations
param_grid <- expand_grid(
  frag_level = names(frag_levels),
  samp_method = samp_methods,
  rep = 1:n_reps
)

pal1 <- c(colorspace::lighten("midnightblue", 0.2), colorspace::lighten("violetred4", 0.2))


# Aggregated state -------------------------------------------------------

# Run all iterations for the aggregated state
aggr_results <- param_grid %>%
  pmap_dfr(~ process_single_rep(
    state = aggr_state,
    frag_level = ..1, 
    samp_method = ..2, 
    rep_id = ..3 
  ))

aggr_results <- readRDS("data-raw/output/aggregated_state_results.rds")

# Averages and confidence intervals for plotting
sSBR_avg <- aggr_results %>%
  unnest(sSBR_main) %>%
  group_by(frag_level, samp_method, effort) %>%
  summarise(
    S_mean = mean(S), 
    S_sd = sd(S),
    n_reps = n(),
    S_se = S_sd / sqrt(n_reps),
    S_low = S_mean - 1.96 * S_se, 
    S_high = S_mean + 1.96 * S_se,
    .groups = "drop"
  ) %>%
  mutate(frag_level = factor(frag_level, levels = c("low", "high")),
         samp_method = factor(samp_method, labels = c("All habitat cells", "Chessboard sampling", "30 random habitat cells")))

dd_avg <- aggr_results %>%
  unnest(dd_main) %>%
  group_by(frag_level, samp_method, distance) %>%
  summarise(
    similarity_mean = mean(similarity),
    similarity_sd = sd(similarity),
    n_reps = n(),
    similarity_se = similarity_sd / sqrt(n_reps),
    similarity_low = similarity_mean - 1.96 * similarity_se,
    similarity_high = similarity_mean + 1.96 * similarity_se,
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
  facet_wrap(~ samp_method, scales = "free_x", ncol = 3) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom")

# Combine: 2 rows x 3 columns = 6 panels
# Single legend at the bottom, shared across all panels
combined_plot <- p1 / p2 + 
  plot_layout(guides = "collect") & 
  guides(color = guide_legend(title = "Fragmentation level", nrow = 1),
         fill = guide_legend(title = "Fragmentation level", nrow = 1)) &
  theme(legend.position = "bottom") #&
  # add additional title on top, without changing the single titles
  # plot_annotation(title = "Single state after ecological dynamics", 
  #                 theme = theme(plot.title = element_text(size = 16, face = "bold")))

print(combined_plot)

ggsave("data-raw/output/aggr_state_combined_plot.png", 
       combined_plot, width = 10, height = 7.5, dpi = 300)



# Random state -----------------------------------------------------------

rand_name <- "State with random species distribution"

rand_results <- param_grid %>%
  pmap_dfr(~ process_single_rep(
    state = rand_state,
    frag_level = ..1, 
    samp_method = ..2, 
    rep_id = ..3 
  ))

saveRDS(rand_results, "data-raw/output/random_state_results.rds")
# rand_results <- readRDS("data-raw/random_state_results.rds")

sSBR_avg_rand <- rand_results %>%
  unnest(sSBR_main) %>%
  group_by(frag_level, samp_method, effort) %>%
  summarise(
    S_mean = mean(S), 
    S_sd = sd(S),
    n_reps = n(),
    S_se = S_sd / sqrt(n_reps),
    S_low = S_mean - 1.96 * S_se, 
    S_high = S_mean + 1.96 * S_se,
    .groups = "drop"
  ) %>%
  mutate(frag_level = factor(frag_level, levels = c("low", "high")),
         samp_method = factor(samp_method, labels = c("All habitat cells", "Chessboard sampling", "30 random habitat cells")))

dd_avg_rand <- rand_results %>%
  unnest(dd_main) %>%
  group_by(frag_level, samp_method, distance) %>%
  summarise(
    similarity_mean = mean(similarity),
    similarity_sd = sd(similarity),
    n_reps = n(),
    similarity_se = similarity_sd / sqrt(n_reps),
    similarity_low = similarity_mean - 1.96 * similarity_se,
    similarity_high = similarity_mean + 1.96 * similarity_se,
    .groups = "drop"
  ) %>%
  mutate(frag_level = factor(frag_level, levels = c("low", "high")),
         samp_method = factor(samp_method, labels = c("All habitat cells", "Chessboard sampling", "30 random habitat cells")))

# Plotting for random state
p1_rand <- ggplot(sSBR_avg_rand, 
             aes(x = effort, y = S_mean, color = frag_level, fill = frag_level)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), alpha = 0.2, color = NA) +
  scale_color_manual(values = pal1) +
  scale_fill_manual(values = colorspace::lighten(pal1, 0.2)) +
  labs(x = "Cumulative area of convex hull", 
       y = "Cumulative species richness",
       title = "Spatially Constrained Rarefaction") +
  facet_wrap(~ samp_method, scales = "free_x", ncol = 3) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom")


p2_rand <- ggplot(dd_avg_rand, 
             aes(x = distance, y = similarity_mean, color = frag_level, fill = frag_level)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = similarity_low, ymax = similarity_high), 
              alpha = 0.2, col = NA) +
  scale_color_manual(values = pal1) +
  scale_fill_manual(values = colorspace::lighten(pal1, 0.2)) +
  labs(x = "Euclidean distance", 
       y = "Similarity (1 - Bray-Curtis dissimilarity)",
       title = "Distance Decay of Similarity") +
  facet_wrap(~ samp_method, scales = "free_x", ncol = 3) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "bottom")

combined_plot_rand <- p1_rand / p2_rand + 
  plot_layout(guides = "collect") & 
  guides(color = guide_legend(title = "Fragmentation level", nrow = 1),
         fill = guide_legend(title = "Fragmentation level", nrow = 1)) &
  theme(legend.position = "bottom"#, 
        # strip.background = element_rect(fill = "slateblue3"),
        # strip.text = element_text(color = "white", size = 12, face = "bold")
        ) &
  plot_annotation(title = "State with random species distribution", 
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

print(combined_plot_rand)

ggsave(paste0("data-raw/output/random_state_combined_plot.png"), 
       combined_plot_rand, width = 15, height = 10, dpi = 300)

# Save not averaged results as csv
saveRDS(rand_results, "data-raw/output/random_state_results.rds")
saveRDS(aggr_results, "data-raw/output/aggregated_state_results.rds")

