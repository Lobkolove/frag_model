source("R/chull_area.R")
source("R/toroidal_clump.R")

data <- data.frame(x = c(1, 3, 3, 1, 2), y = c(1, 1, 3, 3, 2))
chull_area(data, x_col = "x", y_col = "y")

chull_area2(data, x_col = "x", y_col = "y")

picks_hull_area(data, x_col = "x", y_col = "y")

# Test with model output
post_frag_state <- readRDS("output/model_states/0013_ac0.9_frag0.2_hab0.15_r001.rds")$post_fragmentation

# Visualize grid
raster::image(post_frag_state$grid, col = viridis::viridis(100), asp = 1, axes = FALSE, ann = FALSE)
toroidal_clump(post_frag_state$grid)

full_sample <- read.csv("output/sampled_data/0004_ac0.3_frag0.2_hab0.15_r001_samp_all.csv")
chull_area(full_sample)
# chull area much bigger than actual patch size, because of wraparound effect.

chull_area2(data, x_col = "x", y_col = "y")
chull_area2(data, x_col = "x", y_col = "y", torus = TRUE, x_wrap = 4, y_wrap = 4)
chull_area2(full_sample, x_col = "x_loc", y_col = "y_loc", torus = TRUE, x_wrap = 50, y_wrap = 50)


data3 <- data.frame(x = c(0, 1, 4, 5, 0, 1, 4, 5, 0, 1, 4, 5, 0, 1, 4, 5), 
                    y = c(0, 0, 0, 0, 1, 1, 1, 1, 4, 4, 4, 4, 5, 5, 5, 5))

chull_area(data3, x_col = "x", y_col = "y")
toroidal_chull_area(data3, x_col = "x", y_col = "y", Lx = 5, Ly = 5)
torus_hull_area(data3, x_col = "x", y_col = "y", n = 5)
  