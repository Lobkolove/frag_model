# Demo script for distance decay and sSBR functions
library(here)
library(dplyr)
library(vegan)
library(mgcv)
library(scam)

source(here("R", "dist_decay.R"))
source(here("R", "sSBR.R"))
       

# 0. Model sample output --------------------------------------------------

sample_low <- read.csv(here("data-raw/test_Ju_low_rep_1_output_sample.csv"))
sample_high <- read.csv(here("data-raw/test_Ju_high_rep_1_output_sample.csv"))

# At the moment, we're only interested in the samples from right after 
# fragmentation
low41 <- sample_low %>% 
  filter(step == 41)
high41 <- sample_high %>% 
  filter(step == 41)


# 1. Distance decay -------------------------------------------------------

# Compute distance decays
dd_low <- dist_decay(low41,
                     binary = FALSE,
                     method = "bray")
dd_high <- dist_decay(high41,
                      binary = FALSE,
                      method = "bray")

# Plot distance decays
plot(dd_low, col = "darkslategrey")
mtext("Low fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "darkslategrey")
plot(dd_high, col = "steelblue4")
mtext("High fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "steelblue4")


# 3. Spatially constrained sample based rarefaction (sSBR) ----------------

# Compute spatial rarefaction
src_low <- sSBR(low41)
src_high <- sSBR(high41)

# Plot rarefaction curves
plot(src_low, col = "darkslategrey")
mtext("Low fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "darkslategrey")
plot(src_high, col = "steelblue4")
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

