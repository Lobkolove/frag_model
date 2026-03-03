library(dplyr)
library(raster)
source("Model/src/fragmentation.R")
source("Model/src/landscape.R")
source("R/sample_cells.R")
source("R/toroidal_clump.R")
source("R/toroidal_dist.R")
source("R/dist_decay.R")
source("R/sSBR.R")

pre_frag <- readRDS("data-raw/states_frag_0.5_sim_7.rds")$pre_fragmentation

seed_landscape <- pre_frag$master_seed + 1
seed_fragment <- pre_frag$master_seed + 2

low_frag <- fragment(pre_frag, 0.15, 0.25)
high_frag <- fragment(pre_frag, 0.15, 0.75)

low_full <- sample_cells(low_frag, method = "all")
high_full <- sample_cells(high_frag, method = "all")

sSBR_low <- sSBR(low_full, dist_type = "toroidal")
sSBR_high <- sSBR(high_full, dist_type = "toroidal")

# Merge all sSBR objects into a single data frame for plotting with ggplot2
sSBR_low$smooth <- sSBR_low$smooth %>% mutate(fragmentation = "low")
sSBR_high$smooth <- sSBR_high$smooth %>% mutate(fragmentation = "high")
sSBR_smooth <- bind_rows(sSBR_low$smooth, sSBR_high$smooth)

# Plotting
library(ggplot2)
pal2 <- c("pink", "midnightblue")
gg_sSBR_smooth <- ggplot(sSBR_smooth, 
                         aes(x = distance, y = S, color = fragmentation, fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = colorspace::darken(pal2, 0.4)) +
  scale_fill_manual(values = pal2) +
  labs(x = "Toroidal distance between samples",
       y = "Cumulative number of species",
       color = "Level of fragmentation", fill = "Level of fragmentation") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
gg_sSBR_smooth


