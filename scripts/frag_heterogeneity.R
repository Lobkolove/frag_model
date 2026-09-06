# Geometric effect of fragmentation on landscape-level habitat heterogeneity
#
# A single continuous landscape is fragmented at two degrees of fragmentation.
# The top row shows the resulting landscapes, the middle row the set of
# environmental values that survive in each of them, and the bottom row how the
# remaining cells are spread over those values. Habitat amount is held constant,
# so any difference between the two fragmented landscapes is purely geometric: a
# single compact patch samples a narrow part of the environmental gradient but
# holds many cells per value, while many scattered patches sample it more
# broadly and more thinly.

library(raster)
library(viridis)
source("Model/src/landscape.R")
source("Model/src/fragmentation.R")

# Set parameters ---------------------------------------------------------

gr_size <- 25
ac_amount <- 0.75
habitat <- 0.15
fragmentations <- c(0.2, 0.8)

# Number of bins the environmental values are rounded to. Values are rescaled
# to [0, 1] by fbm_fft(), so every cell would otherwise hold a unique value and
# the loss of environmental variety would be invisible. Wider bins also leave
# fewer of them unoccupied in the continuous landscape, whose rescaled extremes
# are isolated single cells.
n_bins <- 50

# Shape of the value swatches in the middle row. NULL gives a square grid
# (8 x 8 for n_bins = 64), a number gives that many columns and as many rows as
# needed (10 with n_bins = 50 gives a flat 10 x 5 block, which keeps the figure
# short). Pick an n_bins that is a multiple of it, or the last row ends short.
swatch_ncol <- 10

# Height of the middle row, relative to the rows of landscapes and histograms
mid_height <- 0.6

# The histograms in the bottom row use coarser bins: a fragmented landscape
# holds too few cells to fill 100 bins, and the resulting comb of single-cell
# spikes hides the difference in bar height the row is meant to show.
n_bins_hist <- 25

# Unique identifier (change for new images)
id <- 2

seed_landscape <- id
seed_fragment <- seed_landscape + 2
# Just to have a different seed for fragmentation


# Landscapes -------------------------------------------------------------

fbm <- fbm_fft(gr_size = gr_size, ac_amount = ac_amount, raster = T, seed = seed_landscape)

# Continuous landscape in the middle, fragmented ones on either side. Both
# fragmented landscapes use the same mask seed, so they are nested versions of
# the same landscape and differ only in how habitat is arranged.
landscapes <- list(
      ls_mask(fbm, habitat = habitat, fragmentation = fragmentations[1], seed = seed_fragment),
      fbm,
      ls_mask(fbm, habitat = habitat, fragmentation = fragmentations[2], seed = seed_fragment)
)


# Colour scale shared by landscapes and value swatches --------------------

# Bin edges and the colour each bin is drawn in. The same mapping is used for
# the rasters and the swatches, so a swatch cell has the colour the matching
# cells have in the landscape above it.
bin_edges <- seq(0, 1, length.out = n_bins + 1)
bin_cols <- magma(n_bins)

# Histogram bins, coloured from the same ramp at their midpoints
hist_edges <- seq(0, 1, length.out = n_bins_hist + 1)
hist_cols <- magma(n_bins_hist)

# Which bins are occupied in a given landscape
occupied_bins <- function(ls) {
      v <- na.omit(raster::values(ls))
      sort(unique(as.integer(cut(v, breaks = bin_edges, include.lowest = TRUE))))
}

# Share of the landscape's cells falling in each histogram bin. Habitat amount
# is held constant, so for the two fragmented landscapes this is cell count up
# to a constant factor: the same cells spread over fewer bins means taller bars.
bin_proportion <- function(ls) {
      v <- na.omit(raster::values(ls))
      as.numeric(table(cut(v, breaks = hist_edges, include.lowest = TRUE))) / length(v)
}

# Draw the occupied bins as a grid of coloured cells. Cells keep the position
# they have in the full palette, so the empty ones make the missing part of the
# environmental gradient visible. Set pack = TRUE to drop the gaps instead.
swatch <- function(bins, n_bins, ncol = NULL, pack = FALSE, empty_col = "grey95") {
      
      n_col <- if (is.null(ncol)) ceiling(sqrt(n_bins)) else ncol
      n_row <- ceiling(n_bins / n_col)
      
      if (pack) {
            cols <- rep(NA, n_bins)
            cols[seq_along(bins)] <- bin_cols[bins]
      } else {
            cols <- rep(NA, n_bins)
            cols[bins] <- bin_cols[bins]
      }
      cols[is.na(cols)] <- empty_col
      
      # Fill the grid left to right, top to bottom
      x <- (seq_len(n_bins) - 1) %% n_col
      y <- n_row - 1 - (seq_len(n_bins) - 1) %/% n_col
      
      plot(NA, xlim = c(0, n_col), ylim = c(0, n_row),
           asp = 1, axes = FALSE, xlab = "", ylab = "")
      rect(x, y, x + 1, y + 1, col = cols, border = "white", lwd = 0.5)
      rect(0, 0, n_col, n_row, border = "black", lwd = 1)
}


# Figure -----------------------------------------------------------------

# Export to file (uncomment to write the figure)
# filename <- paste0("pics/frag_heterogeneity_ac_", ac_amount, "_", id, ".png")
# png(filename, width = 1200, height = 1200)

# Common y scale for the histograms, so bar heights are comparable across panels
props <- lapply(landscapes, bin_proportion)
ymax <- max(unlist(props))

# layout() rather than mfrow, so the middle row can be given less height
layout(matrix(seq_len(9), nrow = 3, byrow = TRUE),
       heights = c(1, mid_height, 1))

par(
      mar   = c(1, 1, 1, 1),   # same for all panels
      oma   = c(2, 2, 2, 2)    # extra outer margins: bottom, left, top, right
)

# Top row: the landscapes
for (i in seq_along(landscapes)) {
      
      raster::image(
            landscapes[[i]],
            asp = 1, axes = FALSE,
            zlim = c(0, 1),   # shared scale, so colours match across panels
            col = bin_cols
      )
      
      # Add masked areas in grey
      raster::image(
            is.na(landscapes[[i]]),
            asp = 1, axes = FALSE,
            col = c(NA, "grey90"),
            add = TRUE
      )
      
      # Add rectangle around the plot
      e <- extent(landscapes[[i]])
      rect(e@xmin, e@ymin, e@xmax, e@ymax, border = "black", lwd = 1)
}

# Middle row: the environmental values still present in each landscape
for (i in seq_along(landscapes)) {
      
      bins <- occupied_bins(landscapes[[i]])
      swatch(bins, n_bins = n_bins, ncol = swatch_ncol)
}

# Bottom row: how the cells are spread over those values. The two fragmented
# landscapes hold the same number of cells, so the narrower distribution is
# also the taller one.
# Wider left margin, to make room for the axis
par(mar = c(1, 4, 1, 1))

for (i in seq_along(landscapes)) {
      
      plot(NA, xlim = c(0, 1), ylim = c(0, ymax),
           axes = FALSE, xlab = "", ylab = "")
      
      rect(hist_edges[-length(hist_edges)], 0, hist_edges[-1], props[[i]],
           col = hist_cols, border = "white", lwd = 0.5)
      
      # Keep the y axis, so the bars are not misread as absolute cell counts
      axis(2, las = 1, cex.axis = 1.1)
      if (i == 1) mtext("proportion of cells", side = 2, line = 2.8, cex = 1.1)
      
      # Add rectangle around the plot
      box(lwd = 1)
}

# dev.off()
