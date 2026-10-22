library(dplyr)
library(ggplot2)
library(raster)
library(viridis)
source("R/sample_cells.R")
source("R/toroidal_clump.R")
source("R/dist_decay.R")
source("R/sSBR.R")

# Load in saved model states for low and high fragmentation
low_frag_state <- readRDS("data-raw/states_frag_0.25_sim_6.rds")$post_fragmentation
high_frag_state <- readRDS("data-raw/states_frag_0.75_sim_6.rds")$post_fragmentation

# Extract landscapes from full states
low_frag_grid <- low_frag_state$grid
high_frag_grid <- high_frag_state$grid

# Visualise landscapes
par(mfrow = c(1, 2), mar = c(1, 2, 1, 1))
image(low_frag_grid, col = viridis(100), asp = 1, axes = FALSE)
image(is.na(low_frag_grid), col = c(NA, "grey90"), asp = 1, add = TRUE)
par(mar = c(1, 1, 1, 2))
image(high_frag_grid, col = viridis(100), asp = 1, axes = FALSE)
image(is.na(high_frag_grid), col = c(NA, "grey90"), asp = 1, add = TRUE)
par(mfrow = c(1, 1))

# Export landscapes to file separately
png("pics/SC_low_frag_grid.png", width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(low_frag_grid, col = viridis(100), asp = 1, axes = FALSE)
image(is.na(low_frag_grid), col = c(NA, "grey90"), asp = 1, add = TRUE)
dev.off()
png("pics/SC_high_frag_grid.png", width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(high_frag_grid, col = viridis(100), asp = 1, axes = FALSE)
image(is.na(high_frag_grid), col = c(NA, "grey90"), asp = 1, add = TRUE)
dev.off()

# Apply sampling methods to low fragmentation state
low_frag_random <- sample_cells(full_state = low_frag_state, method = "random", n_samples = 30)
low_frag_chess <- sample_cells(full_state = low_frag_state, method = "chessboard")
low_frag_full <- sample_cells(full_state = low_frag_state, method = "all")

# Apply sampling methods to high fragmentation state
high_frag_random <- sample_cells(full_state = high_frag_state, method = "random", n_samples = 30)
high_frag_chess <- sample_cells(full_state = high_frag_state, method = "chessboard")
high_frag_full <- sample_cells(full_state = high_frag_state, method = "all")

# Assess alpha and gamma species richness for low and high full samples
# Gamma = ncol of all columns that start with 'sp'
gamma_low <- ncol(low_frag_full %>% select(starts_with("sp")))
gamma_high <- ncol(high_frag_full %>% select(starts_with("sp")))

# Alpha = mean richness across all samples (rows) in each dataset
# for each sample, count species (samples that are not 0) - rowSums not working because we have abundance data
alpha_low <- low_frag_full %>% 
  select(starts_with("sp")) %>% 
  apply(1, function(x) sum(x > 0)) %>% 
  mean()

alpha_high <- high_frag_full %>% 
  select(starts_with("sp")) %>%
  apply(1, function(x) sum(x > 0)) %>%
  mean()

# now alpha for random samples only
alpha_low_random <- low_frag_random %>% 
  select(starts_with("sp")) %>% 
  apply(1, function(x) sum(x > 0)) %>% 
  mean()

alpha_high_random <- high_frag_random %>% 
  select(starts_with("sp")) %>%
  apply(1, function(x) sum(x > 0)) %>%
  mean()

gamma_low_random <- ncol(low_frag_random %>% select(starts_with("sp")))
gamma_high_random <- ncol(high_frag_random %>% select(starts_with("sp")))


# Distance Decay ---------------------------------------------------------

# Compute distance decay for each dataset
dd_low_random <- dist_decay(low_frag_random)
dd_low_chess <- dist_decay(low_frag_chess)
dd_low_full <- dist_decay(low_frag_full)
dd_high_random <- dist_decay(high_frag_random)
dd_high_chess <- dist_decay(high_frag_chess)
dd_high_full <- dist_decay(high_frag_full)

# Merge all distance decay objects into a single data frame for plotting
dd_data_all <- bind_rows(
  dd_low_random$data %>% mutate(fragmentation = "Low", sampling = "Random"),
  dd_low_chess$data %>% mutate(fragmentation = "Low", sampling = "Chessboard"),
  dd_low_full$data %>% mutate(fragmentation = "Low", sampling = "Full"),
  dd_high_random$data %>% mutate(fragmentation = "High", sampling = "Random"),
  dd_high_chess$data %>% mutate(fragmentation = "High", sampling = "Chessboard"),
  dd_high_full$data %>% mutate(fragmentation = "High", sampling = "Full")
) |> 
  mutate(fragmentation = factor(fragmentation, levels = c("Low", "High")),
         sampling = factor(sampling, levels = c("Random", "Chessboard", "Full")))

dd_smooth_all <- bind_rows(
  dd_low_random$smooth %>% mutate(fragmentation = "Low", sampling = "Random"),
  dd_low_chess$smooth %>% mutate(fragmentation = "Low", sampling = "Chessboard"),
  dd_low_full$smooth %>% mutate(fragmentation = "Low", sampling = "Full"),
  dd_high_random$smooth %>% mutate(fragmentation = "High", sampling = "Random"),
  dd_high_chess$smooth %>% mutate(fragmentation = "High", sampling = "Chessboard"),
  dd_high_full$smooth %>% mutate(fragmentation = "High", sampling = "Full")
) |> 
  mutate(fragmentation = factor(fragmentation, levels = c("Low", "High")),
         sampling = factor(sampling, levels = c("Random", "Chessboard", "Full"))) 

# Plot distance decay curves using ggplot2
pal <- c("cyan4", "midnightblue", "violetred4")
gg_dd_full <- ggplot(dd_data_all, aes(x = distance, y = similarity, color = sampling, fill = sampling)) +
    geom_point(alpha = 0.4, shape = 21, color = "#FFFFFF00") +
    geom_line(data = dd_smooth_all, linewidth = 1) +
    geom_ribbon(data = dd_smooth_all, aes(x = distance, ymin = simi_low, ymax = simi_high), 
                alpha = 0.2, inherit.aes = FALSE) +
    scale_color_manual(values = colorspace::darken(pal, 0.4)) +
    scale_fill_manual(values = pal) +
    facet_grid(rows = vars(sampling), 
              cols = vars(fragmentation)) +
    labs(x = "Euclidean distance between samples",
        y = "Similarity (1 - Bray-Curtis dissimilarity)") +
    guides(color = "none", fill = "none") +
    theme_bw()
gg_dd_full

# Export full plot
ggsave("pics/dist_decay_full.png", gg_dd_full, width = 8, height = 6, dpi = 300)

# Plot only averaged curves without points
pal2 <- c("pink", "midnightblue")
gg_dd_smooth <- ggplot(dd_smooth_all, 
                       aes(x = distance, y = similarity, color = fragmentation, fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = simi_low, ymax = simi_high), 
              alpha = 0.4, color = NA) +
  scale_color_manual(values = colorspace::darken(pal2, 0.4)) +
  scale_fill_manual(values = pal2) +
  facet_wrap(~sampling) +
  labs(x = "Euclidean distance between samples",
       y = "Similarity (1 - Bray-Curtis dissimilarity)",
       color = "Level of fragmentation", fill = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "violetred4"),
        strip.text = element_text(color = "white", face = "bold", size = 13),
        legend.position = "bottom")
gg_dd_smooth

# Export smoothed plot
ggsave("pics/dist_decay_smooth.png", gg_dd_smooth, width = 14, height = 8, dpi = 300)


# Spatial Rarefaction Curves ---------------------------------------------

# Compute sSBR for each dataset
sSBR_low_random <- sSBR(low_frag_random)
sSBR_low_chess <- sSBR(low_frag_chess)
sSBR_low_full <- sSBR(low_frag_full)
sSBR_high_random <- sSBR(high_frag_random)
sSBR_high_chess <- sSBR(high_frag_chess)
sSBR_high_full <- sSBR(high_frag_full)

# Merge all sSBR objects into a single data frame for plotting with ggplot2
sSBR_data_all <- bind_rows(
  sSBR_low_random$data %>% mutate(fragmentation = "Low", sampling = "Random"),
  sSBR_low_chess$data %>% mutate(fragmentation = "Low", sampling = "Chessboard"),
  sSBR_low_full$data %>% mutate(fragmentation = "Low", sampling = "Full"),
  sSBR_high_random$data %>% mutate(fragmentation = "High", sampling = "Random"),
  sSBR_high_chess$data %>% mutate(fragmentation = "High", sampling = "Chessboard"),
  sSBR_high_full$data %>% mutate(fragmentation = "High", sampling = "Full")
) |> 
  mutate(fragmentation = factor(fragmentation, levels = c("Low", "High")),
         sampling = factor(sampling, levels = c("Random", "Chessboard", "Full")))

sSBR_smooth_all <- bind_rows(
  sSBR_low_random$smooth %>% mutate(fragmentation = "Low", sampling = "Random"),
  sSBR_low_chess$smooth %>% mutate(fragmentation = "Low", sampling = "Chessboard"),
  sSBR_low_full$smooth %>% mutate(fragmentation = "Low", sampling = "Full"),
  sSBR_high_random$smooth %>% mutate(fragmentation = "High", sampling = "Random"),
  sSBR_high_chess$smooth %>% mutate(fragmentation = "High", sampling = "Chessboard"),
  sSBR_high_full$smooth %>% mutate(fragmentation = "High", sampling = "Full")
) |> 
  mutate(fragmentation = factor(fragmentation, levels = c("Low", "High")),
         sampling = factor(sampling, levels = c("Random", "Chessboard", "Full")))

# Plot sSBR curves using ggplot2
pal <- c("cyan4", "midnightblue", "violetred4")
gg_sSBR_full <- ggplot(sSBR_data_all, 
                       aes(x = distance, y = S, color = sampling)) +
  geom_line(aes(group = id), alpha = 0.1) +
  geom_line(data = sSBR_smooth_all, linewidth = 1) +
  scale_color_manual(values = pal) +
  facet_grid(rows = vars(sampling), 
             cols = vars(fragmentation)) +
  labs(x = "Euclidean distance between samples",
       y = "Cumulative number of species") +
  guides(color = "none", fill = "none") +
  theme_bw()
gg_sSBR_full

# Export full sSBR plot
ggsave("pics/sSBR_full.png", gg_sSBR_full, width = 8, height = 6, dpi = 300)

# Plot only averaged sSBR curves without individual lines
pal2 <- c("pink", "midnightblue")
gg_sSBR_smooth <- ggplot(sSBR_smooth_all, 
                       aes(x = distance, y = S, color = fragmentation, fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = colorspace::darken(pal2, 0.4)) +
  scale_fill_manual(values = pal2) +
  facet_wrap(~sampling) +
  labs(x = "Euclidean distance between samples",
       y = "Cumulative number of species",
       color = "Level of fragmentation", fill = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "violetred4"),
        strip.text = element_text(color = "white", face = "bold", size = 13),
        legend.position = "bottom")
gg_sSBR_smooth

# Export smoothed sSBR plot
ggsave("pics/sSBR_smooth.png", gg_sSBR_smooth, width = 14, height = 8, dpi = 300)

# Toroidal distance ------------------------------------------------------
source("R/toroidal_dist.R")

# Compute distance decay using toroidal distance
dd_low_random_toroidal <- dist_decay(low_frag_random, method = "bray", dist_type = "toroidal")
dd_low_chess_toroidal <- dist_decay(low_frag_chess, method = "bray", dist_type = "toroidal")
dd_low_full_toroidal <- dist_decay(low_frag_full, method = "bray", dist_type = "toroidal")
dd_high_random_toroidal <- dist_decay(high_frag_random, method = "bray", dist_type = "toroidal")
dd_high_chess_toroidal <- dist_decay(high_frag_chess, method = "bray", dist_type = "toroidal")
dd_high_full_toroidal <- dist_decay(high_frag_full, method = "bray", dist_type = "toroidal")

# Merge objects into a single data frame for plotting
dd_toroidal_smooth_all <- bind_rows(
  dd_low_random_toroidal$smooth %>% mutate(fragmentation = "Low", sampling = "Random"),
  dd_low_chess_toroidal$smooth %>% mutate(fragmentation = "Low", sampling = "Chessboard"),
  dd_low_full_toroidal$smooth %>% mutate(fragmentation = "Low", sampling = "Full"),
  dd_high_random_toroidal$smooth %>% mutate(fragmentation = "High", sampling = "Random"),
  dd_high_chess_toroidal$smooth %>% mutate(fragmentation = "High", sampling = "Chessboard"),
  dd_high_full_toroidal$smooth %>% mutate(fragmentation = "High", sampling = "Full")
) |> 
  mutate(fragmentation = factor(fragmentation, levels = c("Low", "High")),
         sampling = factor(sampling, levels = c("Random", "Chessboard", "Full")))

# Plot overall toroidal distance decay curves (low vs high fragmentation) using ggplot2
pal2 <- c("pink", "midnightblue")
gg_dd_toroidal_smooth <- ggplot(dd_toroidal_smooth_all, 
                       aes(x = distance, y = similarity, color = fragmentation, fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = simi_low, ymax = simi_high), 
              alpha = 0.4, color = NA) +
  scale_color_manual(values = colorspace::darken(pal2, 0.4)) +
  scale_fill_manual(values = pal2) +
  facet_wrap(~sampling) +
  labs(x = "Toroidal Distance",
       y = "Similarity (1 - Bray-Curtis Dissimilarity)",
       color = "Level of Fragmentation", fill = "Level of Fragmentation") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "violetred4"),
        strip.text = element_text(color = "white", face = "bold", size = 13),
        legend.position = "bottom")
gg_dd_toroidal_smooth

# Spatial rarefaction with toroidal distance
sSBR_low_random_toroidal <- sSBR(low_frag_random, dist_type = "toroidal")
sSBR_low_chess_toroidal <- sSBR(low_frag_chess, dist_type = "toroidal")
sSBR_low_full_toroidal <- sSBR(low_frag_full, dist_type = "toroidal")
sSBR_high_random_toroidal <- sSBR(high_frag_random, dist_type = "toroidal")
sSBR_high_chess_toroidal <- sSBR(high_frag_chess, dist_type = "toroidal")
sSBR_high_full_toroidal <- sSBR(high_frag_full, dist_type = "toroidal")

# Merge objects into a single data frame for plotting
sSBR_toroidal_smooth_all <- bind_rows(
  sSBR_low_random_toroidal$smooth %>% mutate(fragmentation = "Low", sampling = "Random"),
  sSBR_low_chess_toroidal$smooth %>% mutate(fragmentation = "Low", sampling = "Chessboard"),
  sSBR_low_full_toroidal$smooth %>% mutate(fragmentation = "Low", sampling = "Full"),
  sSBR_high_random_toroidal$smooth %>% mutate(fragmentation = "High", sampling = "Random"),
  sSBR_high_chess_toroidal$smooth %>% mutate(fragmentation = "High", sampling = "Chessboard"),
  sSBR_high_full_toroidal$smooth %>% mutate(fragmentation = "High", sampling = "Full")
) |> 
  mutate(fragmentation = factor(fragmentation, levels = c("Low", "High")),
         sampling = factor(sampling, levels = c("Random", "Chessboard", "Full")))


