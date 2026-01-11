# Demo script for distance decay and sSBR functions
library(here)
library(dplyr)
library(vegan)
library(mgcv)
library(scam)

source(here("R", "dist_decay.R"))
source(here("R", "sSBR.R"))
       

# 0. Model sample output --------------------------------------------------

full_sample <- read.csv(here("data-raw/model_output/fig_3/scale_22_rep_1_output_sample.csv"))

# At the moment, we're only interested in the samples from right after 
# fragmentation
step41 <- model_sample %>% 
  filter(step == 41)

# For the start we'll look at only two levels of fragmentation (for instance 0.2
# and 0.8), to visualise differences between "low" and "high" fragmentation
low_frag <- step41 %>% 
  filter(fragmentation == 0.2)

high_frag <- step41 %>% 
  filter(fragmentation == 0.8)


# 1. Distance decay -------------------------------------------------------

# Compute distance decays
dd_low <- dist_decay(low_frag,
                     binary = FALSE,
                     method = "bray")
dd_high <- dist_decay(high_frag,
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
src_low <- sSBR(low_frag)
src_high <- sSBR(high_frag)

# Plot rarefaction curves
plot(src_low, col = "darkslategrey")
mtext("Low fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "darkslategrey")
plot(src_high, col = "steelblue4")
mtext("High fragmentation", side = 3, line = 0.5, 
      cex = 0.9, col = "steelblue4")
