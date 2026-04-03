library(dplyr)
library(ggplot2)
library(raster)
library(data.table)
library(collapse)
source("Model/parameters.R")
source("Model/src/clean_run.R")
source("Model/src/init_state.R")
source("Model/src/initialize.R")
source("Model/src/landscape.R")
source("Model/src/generate_agents.R")
source("Model/src/distribute_agents.R")
source("Model/src/run_model_step.R")
source("Model/src/birth.R")
source("Model/src/disperse.R")
source("Model/src/death.R")
source("Model/src/immigration.R")
source("Model/src/fragmentation.R")

source("R/sample_cells.R")
source("R/toroidal_clump.R")
source("R/toroidal_dist.R")
source("R/dist_decay.R")
source("R/sSBR.R")


# After ecological dynamics ----------------------------------------------
pre_frag <- readRDS("data-raw/states_frag_0.5_sim_7.rds")$pre_fragmentation

# Apply different fragmentation levels to the same model state
low_frag <- fragment(pre_frag, 0.15, 0.25)
high_frag <- fragment(pre_frag, 0.15, 0.75)

# Full sample for both fragmentation levels
low_full <- sample_cells(low_frag, method = "all")
high_full <- sample_cells(high_frag, method = "all")

pal1 <- c("steelblue4", "magenta3")

  ## sSBR -----------------------------------------------------------------

# Compute convex hull sSBRs for both fragmentation levels
sSBR_area_low <- sSBR(low_full, method = "area")
sSBR_area_high <- sSBR(high_full, method = "area")

# Merge area-based sSBR objects into single data frame for plotting with ggplot2
sSBR_area_low$smooth <- sSBR_area_low$smooth %>% mutate(fragmentation = "low")
sSBR_area_high$smooth <- sSBR_area_high$smooth %>% mutate(fragmentation = "high")
sSBR_area <- bind_rows(sSBR_area_low$smooth, sSBR_area_high$smooth) %>%
  mutate(fragmentation = factor(fragmentation, levels = c("low", "high")))

gg_sSBR_area <- ggplot(sSBR_area, 
                       aes(x = effort, 
                           y = S, 
                           color = fragmentation, 
                           fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = colorspace::lighten(pal2, 0.4)) +
  scale_fill_manual(values = pal2) +
  labs(#title = "Single state after ecological dynamics",
       x = "Cumulative area of convex hull encompassing samples",
       y = "Cumulative number of species",
       color = "Level of fragmentation", fill = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
gg_sSBR_area


# Compute euclidean sSBRs
sSBR_euc_low <- sSBR(low_full, method = "euclidean")
sSBR_euc_high <- sSBR(high_full, method = "euclidean")

# Compute toroidal sSBRs
sSBR_tor_low <- sSBR(low_full, method = "toroidal")
sSBR_tor_high <- sSBR(high_full, method = "toroidal")

# Merge sSBR objects into single data frames for plotting with ggplot2
sSBR_euc_low$smooth <- sSBR_euc_low$smooth %>% mutate(fragmentation = "low")
sSBR_euc_high$smooth <- sSBR_euc_high$smooth %>% mutate(fragmentation = "high")
sSBR_euclidean <- bind_rows(sSBR_euc_low$smooth, sSBR_euc_high$smooth) %>%
  mutate(fragmentation = factor(fragmentation, levels = c("low", "high")))

sSBR_tor_low$smooth <- sSBR_tor_low$smooth %>% mutate(fragmentation = "low")
sSBR_tor_high$smooth <- sSBR_tor_high$smooth %>% mutate(fragmentation = "high")
sSBR_toroidal <- bind_rows(sSBR_tor_low$smooth, sSBR_tor_high$smooth) %>%
  mutate(fragmentation = factor(fragmentation, levels = c("low", "high")))

# Plotting with euclidean distances
gg_sSBR_euc <- ggplot(sSBR_euclidean, 
                      aes(x = effort, 
                          y = S, 
                          color = fragmentation, 
                          fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = colorspace::lighten(pal2, 0.4)) +
  scale_fill_manual(values = pal2) +
  labs(title = "Single state after ecological dynamics (step 40)",
       x = "Euclidean distance between samples",
       y = "Cumulative number of species",
       color = "Level of fragmentation", fill = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
gg_sSBR_euc

# Plotting with toroidal distances
gg_sSBR_tor <- ggplot(sSBR_toroidal, 
                      aes(x = effort, 
                          y = S, 
                          color = fragmentation, 
                          fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = colorspace::lighten(pal2, 0.4)) +
  scale_fill_manual(values = pal2) +
  labs(x = "Toroidal distance between samples",
       y = "Cumulative number of species",
       color = "Level of fragmentation", fill = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
gg_sSBR_tor


  ## distance_decay ------------------------------------------------------

dd_euc_low <- dist_decay(low_full, dist_type = "euclidean")$smooth
dd_euc_high <- dist_decay(high_full, dist_type = "euclidean")$smooth

# merge distance decay objects into single data frame for plotting with ggplot2
dd_euc_low <- dd_euc_low %>% mutate(fragmentation = "low")
dd_euc_high <- dd_euc_high %>% mutate(fragmentation = "high")
dd_euclidean <- bind_rows(dd_euc_low, dd_euc_high) %>%
  mutate(fragmentation = factor(fragmentation, levels = c("low", "high")))

# Plotting distance decay with euclidean distances
gg_dd_euc <- ggplot(dd_euclidean, 
                     aes(x = distance, 
                         y = similarity, 
                         color = fragmentation)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = colorspace::lighten(pal2, 0.4)) +
  labs(#title = "Single state after ecological dynamics",
       x = "Euclidean distance between samples",
       y = "Similarity between samples",
       color = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
gg_dd_euc



# Random Community -------------------------------------------------------

sim_id <- 9
random_state <- clean_run(mod_par = mod_par,
                          var_par = var_par,
                          switch = switch,
                          sim_id = sim_id,
                          record_steps = c("start"),
                          seed = master_seed)

# Export RDS file for random state to use in other scripts
saveRDS(random_state, file = "data-raw/random_state.rds")

random_state_init <- random_state$start

# Apply different fragmentation levels to the random model state
random_low_frag <- fragment(random_state_init, habitat = 0.15, fragmentation = 0.25)
random_high_frag <- fragment(random_state_init, habitat = 0.15, fragmentation = 0.75)

# Sample all cells for both fragmentation levels
random_low_full <- sample_cells(random_low_frag, method = "all")
random_high_full <- sample_cells(random_high_frag, method = "all")

pal2 <- c("midnightblue", "violetred3")

  ## Convex hull area ####

# Compute area-based sSBRs for random state
random_sSBR_area_low <- sSBR(random_low_full, method = "area", cutoff = FALSE)
random_sSBR_area_high <- sSBR(random_high_full, method = "area", cutoff = FALSE)

# Plotting area-based sSBR for random state
random_sSBR_area_low$smooth <- random_sSBR_area_low$smooth %>% mutate(fragmentation = "low")
random_sSBR_area_high$smooth <- random_sSBR_area_high$smooth %>% mutate(fragmentation = "high")
random_sSBR_area <- bind_rows(random_sSBR_area_low$smooth, random_sSBR_area_high$smooth) %>%
  mutate(fragmentation = factor(fragmentation, levels = c("low", "high")))

gg_random_sSBR_area <- ggplot(random_sSBR_area, 
                              aes(x = effort, 
                                  y = S, 
                                  color = fragmentation, 
                                  fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = colorspace::lighten(pal2, 0.4)) +
  scale_fill_manual(values = pal2) +
  labs(#title = "Random distribution of species across the landscape",
       x = "Cumulative area of convex hull encompassing samples",
       y = "Cumulative number of species",
       color = "Level of fragmentation", fill = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
gg_random_sSBR_area

  ## Distance decay ####

random_dd_euc_low <- dist_decay(random_low_full, dist_type = "euclidean")$smooth
random_dd_euc_high <- dist_decay(random_high_full, dist_type = "euclidean")$smooth

# Merge distance decay objects into single data frame for plotting with ggplot2
random_dd_euc_low <- random_dd_euc_low %>% mutate(fragmentation = "low")
random_dd_euc_high <- random_dd_euc_high %>% mutate(fragmentation = "high")
random_dd_euclidean <- bind_rows(random_dd_euc_low, random_dd_euc_high) %>%
  mutate(fragmentation = factor(fragmentation, levels = c("low", "high")))

# Plotting distance decay with euclidean distances for random state
gg_random_dd_euc <- ggplot(random_dd_euclidean, 
                           aes(x = distance, 
                               y = similarity, 
                               color = fragmentation)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = colorspace::lighten(pal2, 0.4)) +
  labs(#title = "Random distribution of species across the landscape",
       x = "Euclidean distance between samples",
       y = "Similarity between samples",
       color = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
gg_random_dd_euc


  ## Distance ####

# Compute euclidean sSBRs for random state
random_sSBR_euc_low <- sSBR(random_low_full, cutoff = FALSE)
random_sSBR_euc_high <- sSBR(random_high_full, cutoff = FALSE)

# Compute toroidal sSBRs for random state
random_sSBR_tor_low <- sSBR(random_low_full, method = "toroidal")
random_sSBR_tor_high <- sSBR(random_high_full, method = "toroidal")

# Plotting sSBR for random state
random_sSBR_euc_low$smooth <- random_sSBR_euc_low$smooth %>% mutate(fragmentation = "low")
random_sSBR_euc_high$smooth <- random_sSBR_euc_high$smooth %>% mutate(fragmentation = "high")
random_sSBR_euclidean <- bind_rows(random_sSBR_euc_low$smooth, random_sSBR_euc_high$smooth) %>%
  mutate(fragmentation = factor(fragmentation, levels = c("low", "high")))

random_sSBR_tor_low$smooth <- random_sSBR_tor_low$smooth %>% mutate(fragmentation = "low")
random_sSBR_tor_high$smooth <- random_sSBR_tor_high$smooth %>% mutate(fragmentation = "high")
random_sSBR_toroidal <- bind_rows(random_sSBR_tor_low$smooth, random_sSBR_tor_high$smooth) %>%
  mutate(fragmentation = factor(fragmentation, levels = c("low", "high")))

gg_random_sSBR_euc <- ggplot(random_sSBR_euclidean, 
                              aes(x = effort, 
                                  y = S, 
                                  color = fragmentation, 
                                  fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = colorspace::lighten(pal2, 0.4)) +
  scale_fill_manual(values = pal2) +
  labs(title = "Single state with random species distribution (step 0)",
       x = "Euclidean distance between samples",
       y = "Cumulative number of species",
       color = "Level of fragmentation", fill = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
gg_random_sSBR_euc

gg_random_sSBR_tor <- ggplot(random_sSBR_toroidal, 
                              aes(x = effort, 
                                  y = S, 
                                  color = fragmentation, 
                                  fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = colorspace::lighten(pal2, 0.4)) +
  scale_fill_manual(values = pal2) +
  labs(title = "Random distribution of species across the landscape",
       x = "Toroidal distance between samples",
       y = "Cumulative number of species",
       color = "Level of fragmentation", fill = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
gg_random_sSBR_tor



# Plots to export --------------------------------------------------------

# For single state after ecological dynamics, 
# both convex hull and distance-decay, side by side
# Single title for both plots, legend at the bottom (suppress individual plot titles)
library(patchwork)
gg_sSBR_area <- gg_sSBR_area + guides(fill = "none")
gg_sSBR_area + gg_dd_euc + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Single state after ecological dynamics (step 40)",
                  theme = theme(plot.title = element_text(hjust = 0.1, size = 16))) & 
  theme(legend.position = "bottom", plot.title.position = "panel")
ggsave("pics/sim_state.png", width = 12, height = 6, dpi = 300)

# For random community, both convex hull and distance-decay, side by side
gg_random_sSBR_area <- gg_random_sSBR_area + guides(fill = "none")
gg_random_sSBR_area + gg_random_dd_euc +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Single state with random species distribution (step 0)",
                  theme = theme(plot.title = element_text(hjust = 0.1, size = 16))) & 
  theme(legend.position = "bottom", plot.title.position = "panel")
ggsave("pics/random_state.png", width = 12, height = 6, dpi = 300)

# Euclidean sSBRs for comparison 
gg_sSBR_euc
ggsave("pics/sim_state_euc.png", width = 12, height = 6, dpi = 300)
gg_random_sSBR_euc
ggsave("pics/random_state_euc.png", width = 12, height = 6, dpi = 300)

# without Cutoff
gg_random_sSBR_area + gg_random_sSBR_euc +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Random distribution of species across the landscape (no cutoff)",
                  theme = theme(plot.title = element_text(hjust = 0.1, size = 16))) & 
  theme(legend.position = "bottom", plot.title.position = "panel")
ggsave("pics/random_state_no_cutoff.png", width = 12, height = 6, dpi = 300)
