library(raster)
library(viridis)
library(ggplot2)
source("Model/src/landscape.R")
source("Model/src/fragmentation.R")
source("R/tile_raster.R")

# Set default parameters
gr_size <- 50
ac_amount <- 0.5

# Unique identifier (change for new images)
id <- 81

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

fragmentation <- 0.5
seed_fragment <- seed_landscape + 4
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


# 
# Continuous vs. fragmented landscapes (2x2) -----------------------------

# Same underlying landscape shown continuous (A) and at three levels of
# fragmentation (B-D), which trade a single large patch for many small ones.

annotations <- c("A", "B", "C", "D")
fragmentations <- c(0.1, 0.4, 0.7)  # B, C, D

ac_amount <- 0.7
seed_landscape <- 20
seed_fragment <- seed_landscape + 1

fbm <- fbm_fft(gr_size = gr_size, ac_amount = ac_amount, raster = T, seed = seed_landscape)

# Export to file on Desktop
# png(filename, width = 800, height = 800)
# par(
#       mfrow = c(2, 2),
#       mar   = c(2, 2, 2, 2),   # same for all panels, top margin leaves room for labels
#       oma   = c(1, 1, 1, 1)    # extra outer margins: bottom, left, top, right
# )

for (i in seq_along(annotations)) {
      
      # Panel A shows the continuous landscape, B-D the fragmented versions
      if (i == 1) {
            panel_ls <- fbm
      } else {
            panel_ls <- ls_mask(fbm, habitat = 0.15,
                                fragmentation = fragmentations[i - 1],
                                seed = seed_fragment)
      }
      
      filename <- paste0("pics/loss_frag_2x2_", annotations[i], ".png")
      png(filename, width = 600, height = 600)
      par(mar = c(0,0,0,0))
  
      # Draw raster
      raster::image(
            panel_ls,
            asp = 1, axes = FALSE,
            col = viridis(100)
      )
      
      # Add masked areas in grey
      if (i > 1) {
            raster::image(
                  is.na(panel_ls),
                  asp = 1, axes = FALSE,
                  xaxt = "n", yaxt = "n",
                  col = c(NA, "grey90"),
                  add = TRUE
            )
      }
      
      # Add rectangle around the plot
      e <- extent(panel_ls)
      rect(e@xmin, e@ymin, e@xmax, e@ymax, border = "black", lwd = 1)
      
      dev.off()
}

dev.off()


# Distribution of environmental values -----------------------------------

# Average histogram of cell values across many replicate landscapes, for the
# same levels of autocorrelation as in the multiple-landscapes figure.
# Values are rescaled to [0, 1] within each landscape by fbm_fft().

annotations <- c("A", "B", "C")
ac_amounts <- c(1, 0.5, 0)
n_reps <- 50

breaks <- seq(0, 1, length.out = 41)
mids <- breaks[-length(breaks)] + diff(breaks) / 2

# One replicate landscape -> density per bin
rep_density <- function(ac, seed) {
      ls_rep <- fbm_fft(gr_size = gr_size, ac_amount = ac, raster = FALSE, seed = seed)
      hist(ls_rep, breaks = breaks, plot = FALSE)$density
}

# Replicates share seeds across ac levels, so each level sees the same set of
# noise realisations and only the spectral filtering differs
dens_mean <- matrix(NA_real_, nrow = length(mids), ncol = length(ac_amounts))
dens_sd   <- matrix(NA_real_, nrow = length(mids), ncol = length(ac_amounts))

for (j in seq_along(ac_amounts)) {
      reps <- sapply(seq_len(n_reps), function(r) rep_density(ac_amounts[j], seed = id * 1000 + r))
      dens_mean[, j] <- rowMeans(reps)
      dens_sd[, j]   <- apply(reps, 1, sd)
}

ymax <- max(dens_mean + dens_sd)

# Export to file (uncomment to write the figure)
filename <- paste0("pics/fbm_value_distributions_", id, ".png")
png(filename, width = 1200, height = 800)

par(
      mfrow = c(1, 3),
      mar   = c(4, 4, 3, 1),   # same for all panels
      oma   = c(0.1, 2, 2, 2)    # extra outer margins: bottom, left, top, right
)

for (j in seq_along(ac_amounts)) {
      
      plot(NA, xlim = c(0, 1), ylim = c(0, ymax),
           xlab = "Environmental value", ylab = "Density",
           las = 1, cex.lab = 1.4, cex.axis = 1.2)
      
      # Mean histogram across replicates
      rect(breaks[-length(breaks)], 0, breaks[-1], dens_mean[, j],
           col = viridis(length(ac_amounts))[j], border = "white")
      
      # Variability across replicates
      segments(mids, pmax(dens_mean[, j] - dens_sd[, j], 0),
               mids, dens_mean[, j] + dens_sd[, j],
               col = "grey30", lwd = 1)
      
      # Add bold subplot label above the plot, left-aligned
      mtext(
            annotations[j],
            side = 3, line = 1, adj = 0.05, cex = 2, font = 2
      )
      
      # Report the level of autocorrelation shown
      mtext(
            bquote(ac == .(ac_amounts[j])),
            side = 3, line = 1, adj = 0.95, cex = 1.2
      )
}

dev.off()


# Animated landscape series ----------------------------------------------

# Same idea as the GIF in docs/early_extensions.qmd, on a 50 x 50 grid: one
# animation sweeping the autocorrelation range, one sweeping the fragmentation
# range. Both keep their seed fixed, so consecutive frames are the same
# underlying random field seen at a different setting rather than a new
# landscape each time, and the series morphs smoothly.

library(gifski)   # install.packages("gifski") if not available

id <- 20
gr_size <- 50
ac_amount <- 0.7
habitat <- 0.15
seed_landscape <- 81
seed_fragment <- 21

# Number of frames per sweep. Finer steps give a smoother but heavier GIF.
step <- 0.1

# Run the sweep back down again, so the animation loops without a jump
ping_pong <- function(x) c(x, rev(x)[-c(1, length(x))])

# One frame, drawn in the same style as the static figures above
draw_frame <- function(ls, label) {
      
      par(mar = c(0, 0, 3, 0))
      
      image(ls, asp = 1, axes = FALSE,
            zlim = c(0, 1),   # fixed scale, so colours do not shift between frames
            col = viridis(100))
      
      # Add masked areas in grey
      if (any(is.na(raster::values(ls)))) {
            image(is.na(ls), asp = 1, axes = FALSE,
                  col = c(NA, "grey90"), add = TRUE)
      }
      
      # Add rectangle around the plot
      e <- extent(ls)
      rect(e@xmin, e@ymin, e@xmax, e@ymax, border = "black", lwd = 1)
      
      mtext(label, side = 3, line = 0.5, cex = 4)
}

      ## Varying autocorrelation -----

ac_levels <- ping_pong(seq(0, 1, by = step))

filename <- paste0("pics/fbm_autocorrelation_", id, ".gif")
save_gif(
      for (ac in ac_levels) {
            ls_ac <- fbm_fft(gr_size = gr_size, ac_amount = ac,
                             raster = T, seed = seed_landscape)
            draw_frame(ls_ac, sprintf("autocorrelation = %.1f", ac))
      },
      gif_file = filename,
      width = 600, height = 660, delay = 0.5, progress = FALSE
)

      ## Varying fragmentation -----

frag_levels <- ping_pong(seq(0, 1, by = step))

fbm <- fbm_fft(gr_size = gr_size, ac_amount = ac_amount, raster = T, seed = seed_landscape)

filename <- paste0("pics/fbm_fragmentation_", id, ".gif")
save_gif(
      for (fr in frag_levels) {
            ls_fr <- ls_mask(fbm, habitat = habitat,
                             fragmentation = fr, seed = seed_fragment)
            draw_frame(ls_fr, sprintf("fragmentation = %.1f", fr))
      },
      gif_file = filename,
      width = 600, height = 660, delay = 0.5, progress = FALSE
)
