################################################################################
#
# Establishment probability under short vs random dispersal on continuous
# landscapes
#
################################################################################
#
# Average probability that an offspring passes the niche filter at its
# destination, for short dispersal at mean distances 1, 2, 4 and 8 and for
# random habitat dispersal, on continuous landscapes with autocorrelation
# 0, 0.5 and 1.
#
# No populations are simulated. Each trial is one random individual: a random
# species optimum `u`, a random start cell, one dispersal step, and the model's
# survival probability exp(-(e - u)^2 / (2 nb^2)) at the destination. On a
# continuous landscape every cell is habitat and there is no crowding, so this
# is the whole establishment cascade apart from the constant birth rate.
#
# The parent is first passed through the same niche filter at its start cell.
# Without that step, the parent's cell would be unrelated to its optimum and
# every dispersal mode would give the same average by construction; with it,
# the parent is a plausible established individual and short dispersal can
# benefit from landing in an environment correlated with the parent's.
# `CONDITION_ON_PARENT <- FALSE` gives the unconditioned baseline.
#
################################################################################

library(here)
library(data.table)
library(checkmate)
library(withr)
library(scales)

source(here("Model/src/landscape.R"))


# Configuration ---------------------------------------------------------------

AC_LEVELS <- c(0, 0.5, 1)
DISP_DISTS <- c(1, 2, 4, 8)

GRID_SIZE <- 50
N_SPECIES <- 1000
NICHE_BREADTH <- 0.1

N_LANDSCAPES <- 10L   # landscapes per ac level
N_TRIALS <- 50000L    # individuals per landscape and dispersal mode
CONDITION_ON_PARENT <- TRUE
SEED <- 42L


# Dispersal -------------------------------------------------------------------

# Vectorised equivalents of the model's `toroidal_disperse()` (exponential
# kernel) and `random_disperse(force_habitat = TRUE)` on a continuous grid.
# Positions are (row, col) matrices.
disperse_short <- function(loc, d_mean, grid_size) {
  n <- nrow(loc)
  distance <- rexp(n, 1 / d_mean)
  direction <- runif(n, 0, 2 * pi)
  raw <- loc + cbind(round(cos(direction) * distance), round(sin(direction) * distance))
  ((raw - 1) %% grid_size) + 1
}

disperse_random <- function(n, grid_size) {
  cbind(sample.int(grid_size, n, replace = TRUE), sample.int(grid_size, n, replace = TRUE))
}


# Probe -----------------------------------------------------------------------

#' Mean niche-survival probability at the destination for one landscape.
#'
#' @return A `data.table` with one row per dispersal mode.
probe_landscape <- function(env, n_trials, nb, n_species, condition_on_parent) {

  grid_size <- nrow(env)
  survival <- function(e, u) exp(-(e - u)^2 / (2 * nb^2))

  u <- seq(0, 1, length.out = n_species)[sample.int(n_species, n_trials, replace = TRUE)]
  start <- disperse_random(n_trials, grid_size)

  if (condition_on_parent) {
    keep <- survival(env[start], u) > runif(n_trials)
    u <- u[keep]
    start <- start[keep, , drop = FALSE]
  }
  n <- length(u)

  modes <- c(paste0("short_d", DISP_DISTS), "random")
  rbindlist(lapply(modes, function(m) {
    dest <- if (m == "random") {
      disperse_random(n, grid_size)
    } else {
      disperse_short(start, as.numeric(sub("short_d", "", m)), grid_size)
    }
    data.table(
      mode = m,
      n = n,
      p_establish = mean(survival(env[dest], u))
    )
  }))
}


# Analysis --------------------------------------------------------------------

set.seed(SEED)

results <- rbindlist(lapply(AC_LEVELS, function(ac) {
  rbindlist(lapply(seq_len(N_LANDSCAPES), function(r) {
    env <- fbm_fft(gr_size = GRID_SIZE, ac_amount = ac, raster = FALSE)
    cbind(
      ac = ac,
      landscape = r,
      probe_landscape(env, N_TRIALS, NICHE_BREADTH, N_SPECIES, CONDITION_ON_PARENT)
    )
  }))
}))

results[, mode := factor(mode, levels = c(paste0("short_d", DISP_DISTS), "random"))]

summary_tbl <- dcast(
  results[, .(p = mean(p_establish)), by = .(ac, mode)],
  ac ~ mode,
  value.var = "p"
)

cat("\nMean establishment probability (niche pass at destination), averaged over",
    N_LANDSCAPES, "landscapes per ac level\n",
    "parents conditioned on passing the niche filter at their start cell:",
    CONDITION_ON_PARENT, "\n\n")
print(summary_tbl, digits = 3)

fwrite(results, here("output/establishment_probs_continuous.csv"))
cat("\nPer-landscape results written to output/establishment_probs_continuous.csv\n")
