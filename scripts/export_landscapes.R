library(raster)
library(viridis)
library(ggplot2)
source("Model/src/landscape.R")
source("Model/src/fragmentation.R")
source("R/tile_raster.R")

# Set parameters
gr_size <- 50
ac_amount <- 0.7

# Unique identifier (change for new images)
id <- 151

# Single full landscape ---------------------------------------------------------

seed_landscape <- id

fbm <- fbm_fft(gr_size = gr_size, ac_amount = ac_amount, raster = T, seed = seed_landscape)
# As raster to have consistent plotting with ls_mask output

# Plot to screen
par(mar = c(0, 0, 0, 0))
image(fbm, asp = 1, axes = FALSE,
      col = viridis(100))
dev.off()

# Export to file
filename <- paste0("pics/fbm_full_ac_", ac_amount, "_", id, ".png")
png(filename, width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(fbm, asp = 1, axes = FALSE,
      col = viridis(100))
dev.off()


# Tiled landscape --------------------------------------------------------

tiled <- tile_raster(fbm, n = 3)

# Plot to screen
par(mar = c(0, 0, 0, 0))
image(tiled, asp = 1, axes = FALSE,
      col = viridis(100))
dev.off()

# Export to file
filename <- paste0("pics/fbm_tiled_ac_", ac_amount, "_", id, ".png")
png(filename, width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(tiled, asp = 1, axes = FALSE,
      col = viridis(100))
dev.off()


# Fragmented landscape ---------------------------------------------------

fragmentation <- 0.2
seed_fragment <- seed_landscape + 1 
# Just to have a different seed for fragmentation

frag <- ls_mask(fbm, habitat = 0.15, fragmentation = fragmentation, seed = seed_fragment)

# Plot to screen
par(mar = c(0, 0, 0, 0))
image(frag, asp = 1, axes = FALSE,
      col = viridis(100))
# Add masked areas in grey
image(is.na(frag), asp = 1, axes = FALSE,
      col = c(NA, "grey90"), add = TRUE)
dev.off()

# Export to file
filename <- paste0("pics/fbm_fragmented_ac_", ac_amount, "_frag_", fragmentation, "_", id, ".png")
png(filename, width = 1200, height = 1200)
par(mar = c(0, 0, 0, 0))
image(frag, asp = 1, axes = FALSE,
      col = viridis(100))
# Add masked areas in grey
image(is.na(frag), asp = 1, axes = FALSE,
      col = c(NA, "grey90"), add = TRUE)
dev.off()


# Multiple landscapes ----------------------------------------------------

annotations <- c("A", "B", "C")

      ## Varying autocorrelation -----
ac_amounts <- c(1, 0.5, 0)
seed_landscape <- id

filename <- paste0("pics/fbm_multiple_ac_", id, ".png")
png(filename, width = 1200, height = 400)
par(mfrow = c(1, 3), 
      mar   = c(1, 1, 1, 1),   # same for all panels
      oma   = c(0.1, 2, 4, 2)    # extra outer margins: bottom, left, top, right
)
for (ac in ac_amounts) {
      fbm <- fbm_fft(gr_size = gr_size, ac_amount = ac, raster = T, seed = seed_landscape)
      image(fbm, asp = 1, axes = FALSE,
            col = viridis(100))
      
      # Add rectangle around the plot
      e <- extent(fbm)
      rect(e@xmin, e@ymin, e@xmax, e@ymax, border = "black", lwd = 1)

      # Add bold subplot label above the plot, left-aligned
      mtext(
            annotations[which(ac_amounts == ac)],
            side = 3, line = 1, adj = 0.05, cex = 2, font = 2
      )
}
dev.off()

      ## Varying fragmentation -----

fragmentations <- c(0.2, 0.5, 0.8)
seed_fragment <- 21

ac_amount <- 0.7
fbm <- fbm_fft(gr_size = gr_size, ac_amount = ac_amount, raster = T, seed = seed_landscape)

filename <- paste0("pics/fbm_multiple_frag_", id, ".png")

png(filename, width = 1200, height = 400)
par(
      mfrow = c(1, 3),
      mar   = c(1, 1, 1, 1),   # same for all panels
      oma   = c(0.1, 2, 4, 2)    # extra outer margins: bottom, left, top, right
)
for (frag in fragmentations) {
      
      frag_ls <- ls_mask(fbm, habitat = 0.15, fragmentation = frag, seed = seed_fragment)
      
      # Draw raster
      raster::image(
            frag_ls,
            asp = 1, axes = FALSE,
            col = viridis(100)
      )

      # Add masked areas in grey
      raster::image(
            is.na(frag_ls),
            asp = 1, axes = FALSE,
            col = c(NA, "grey90"),
            add = TRUE
      )

      # Add rectangle around the plot
      e <- extent(frag_ls)
      rect(e@xmin, e@ymin, e@xmax, e@ymax, border = "black", lwd = 1)

      # Add bold subplot label above the plot, left-aligned 
      mtext(
            annotations[which(fragmentations == frag)],
            side = 3, line = 1, adj = 0.05, cex = 2, font = 2
      )
}
dev.off()
      