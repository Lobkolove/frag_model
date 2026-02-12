# Demo script for distance decay and sSBR functions
library(dplyr)
library(vegan)
library(mgcv)
library(scam)

source("R/sample_cells.R")
source("R/toroidal_clump.R")
source("R/dist_decay.R")
source("R/sSBR.R")
source("R/toroidal_dist.R")


# 0. Setup ----------------------------------------------------------------

# Read in full states object and extract post-fragmentation states
states_high <- readRDS("data-raw/states_frag_0.75_sim_6.rds")
states_low <- readRDS("data-raw/states_frag_0.25_sim_6.rds")

post_frag_high <- states_high$post_fragmentation
post_frag_low <- states_low$post_fragmentation

# Create samples for distance decay and sSBR
high_rand <- sample_cells(post_frag_high, method = "random", n_samples = 30)
low_rand <- sample_cells(post_frag_low, method = "random", n_samples = 30)

# 1. Distance decay -------------------------------------------------------
  
# Compute distance decays
dd_low_eu <- dist_decay(low_rand, method = "bray")
dd_low_tor <- dist_decay(low_rand, method = "bray", dist_type = "toroidal")
dd_high_eu <- dist_decay(high_rand, method = "bray")
dd_high_tor <- dist_decay(high_rand, method = "bray", dist_type = "toroidal")

# Plot distance decays
par(mfrow = c(2, 2))
plot(dd_low_eu, col = "darkslategrey")
mtext("Low fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "darkslategrey")
plot(dd_high_eu, col = "steelblue4")
mtext("High fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "steelblue4")
plot(dd_low_tor, col = "darkslategrey")
mtext("Low fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "darkslategrey")
plot(dd_high_tor, col = "steelblue4")
mtext("High fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "steelblue4")  


# 2. Spatially constrained sample based rarefaction (sSBR) ----------------

# Compute spatial rarefaction
src_low_eu <- sSBR(low_rand, dist_type = "euclidean")
src_low_tor <- sSBR(low_rand, dist_type = "toroidal")
src_high_eu <- sSBR(high_rand, dist_type = "euclidean")
src_high_tor <- sSBR(high_rand, dist_type = "toroidal")

# Plot rarefaction curves
par(mfrow = c(2, 2))
plot(src_low_eu, col = "darkslategrey")
mtext("Low fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "darkslategrey")
plot(src_high_eu, col = "steelblue4")
mtext("High fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "steelblue4")
plot(src_low_tor, col = "darkslategrey")
mtext("Low fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "darkslategrey")
plot(src_high_tor, col = "steelblue4")
mtext("High fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "steelblue4")


# Example with old sample ----
full_sample <- read.csv(here("data-raw/model_output/fig_3/scale_22_rep_1_output_sample.csv")) %>% 
  filter(step == 41)

f_low <- full_sample %>% 
  filter(fragmentation == 0.2)
f_high <- full_sample %>% 
  filter(fragmentation == 0.8)

dd_low <- dist_decay(f_low,
                     binary = FALSE,
                     method = "bray")
dd_high <- dist_decay(f_high,
                      binary = FALSE,
                      method = "bray")

# Plot distance decays
plot(dd_low, col = "darkslategrey")
mtext("Low fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "darkslategrey")
plot(dd_high, col = "steelblue4")
mtext("High fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "steelblue4")

src_low <- sSBR(f_low)
src_high <- sSBR(f_high)

plot(src_low, col = "darkslategrey")
mtext("Low fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "darkslategrey")
plot(src_high, col = "steelblue4")
mtext("High fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "steelblue4")


# Example with different sampling methods ---------------------------------



# Compute distance decays
dd_ran <- dist_decay(sample_ran, binary = FALSE, method = "bray")
dd_full <- dist_decay(sample_full, binary = FALSE, method = "bray")
dd_cb <- dist_decay(sample_cb, binary = FALSE, method = "bray")

# Plot distance decays
plot(dd_ran, col = "darkslategrey")
mtext("30 Random Samples", side = 3, line = 0.5, cex = 0.9, col = "darkslategrey")
plot(dd_full, col = "steelblue4")
mtext("Full grid", side = 3, line = 0.5, cex = 0.9, col = "steelblue4")
plot(dd_cb, col = "violetred4")
mtext("Chessboard", side = 3, line = 0.5, cex = 0.9, col = "violetred4")


# Compute SRCs
src_ran <- sSBR(sample_ran)
src_full <- sSBR(sample_full)
src_cb <- sSBR(sample_cb)

# Plot distance decays
plot(src_ran, col = "darkslategrey")
mtext("30 Random Samples", side = 3, line = 0.5, cex = 0.9, col = "darkslategrey")
plot(src_full, col = "steelblue4")
mtext("Full grid", side = 3, line = 0.5, cex = 0.9, col = "steelblue4")
plot(src_cb, col = "violetred4")
mtext("Chessboard sampling", side = 3, line = 0.5, cex = 0.9, col = "violetred4")


# Full iterative ----------------------------------------------------------

# Read in data
sampling_methods <- c("ran", "full", "cb")
fragmentation_levels <- c("low", "high")

data <- list()

for (frag in fragmentation_levels) {
  for (samp in sampling_methods) {
    key <- paste(frag, samp, sep = "_")
    file <- paste0("data-raw/post_frag_", frag, "_sample_", samp, ".csv")
    data[[key]] <- read.csv(file)
  }
}

# Compute distance decays
dd <- list()
for (samp in sampling_methods) {
  dd[[paste0("low_", samp)]]  <- dist_decay(data[[paste0("low_", samp)]],  method = "bray")
  dd[[paste0("high_", samp)]] <- dist_decay(data[[paste0("high_", samp)]], method = "bray")
}

# Plot curves
cols <- c(low = "darkslategrey", high = "steelblue4")

for (samp in sampling_methods) {
  
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  
  # low fragmentation
  plot(
    dd[[paste0("low_", samp)]],
    col = cols["low"]
  )
  mtext("Low fragmentation", side = 3, line = 0.3, cex = 0.9)
  mtext(
    paste("Sampling method:", samp),
    side = 1, line = 3, cex = 0.8, col = "grey30"
  )
  
  # high fragmentation
  plot(
    dd[[paste0("high_", samp)]],
    col = cols["high"]
  )
  mtext("High fragmentation", side = 3, line = 0.3, cex = 0.9)
  mtext(
    paste("Sampling method:", samp),
    side = 1, line = 3, cex = 0.8, col = "grey30"
  )
}


