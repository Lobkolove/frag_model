################################################################################
#
# Establishment success probabilities under different dispersal modes
#
################################################################################
#
# Purpose
# -------
# Decompose the per-capita probability that a birth attempt results in an
# established offspring, separately for each dispersal mode, and attribute the
# failures to the stage of the cascade that caused them.
#
# Can establishment be probed without running the whole model?
# ------------------------------------------------------------
# Yes, exactly. Within a single call to `birth()` the establishment cascade
# reads only three things:
#
#   1. `grid`        - to test habitat vs matrix and to read the environmental
#                      value `e` at the destination;
#   2. `agents`      - to count occupants of the destination cell for the
#                      carrying-capacity checks;
#   3. the parent    - its `species_id` (-> `u`, `nb`) and its location.
#
# `birth()` never mutates the grid, and the only within-step feedback is that
# newborns added earlier in the loop count towards the capacity of later
# attempts. So conditional on (grid, agents, species_par), every birth attempt
# is independent, and replaying attempts against a frozen state reproduces the
# real per-attempt probabilities exactly. (The newborn feedback is a
# second-order effect; `crowding_headroom()` below quantifies how much slack
# there is, and it turns out to be irrelevant at these densities.)
#
# What CANNOT be obtained standalone is the *composition* of `agents` - which
# species are present and how well matched they are to their cells. That
# composition is itself the main thing the dispersal mode changes. So the
# script deliberately crosses both probes over both states:
#
#            probe as "random"   probe as "short_long"
#   state A       (a)                    (b)          <- from a short run
#   state B       (c)                    (d)          <- from a random run
#
# Comparing (a) vs (b) isolates the *mode* effect holding composition fixed.
# Comparing (a) vs (c) isolates the *composition* effect holding mode fixed.
#
# The probe reuses the model's own `random_disperse()` and `toroidal_disperse()`
# rather than reimplementing them, so it validates the real dispersal chain.
#
################################################################################

library(here)
library(data.table)
library(raster)
library(collapse)

source(here("Model/src/disperse.R"))


# Helpers ---------------------------------------------------------------------

#' Reconstruct species parameters from a recorded state's metadata.
#'
#' Under `switch$species_specific_par == 0` (the setting all logged runs used)
#' `species_par` is fully deterministic: `n_value` is an evenly spaced sequence
#' on [0, 1] and every other parameter is constant. Rebuilding it from `meta`
#' rather than sourcing `parameters.R` keeps the probe independent of whatever
#' the current working parameter file happens to say.
species_par_from_meta <- function(
  meta,
  n_species = 1000,
  birth_rate = 0.85,
  death_rate = 0.25
) {
  data.table(
    species_id = 1:n_species,
    n_value = seq(from = 0, to = 1, length.out = n_species),
    birth_rate = birth_rate,
    death_rate = death_rate,
    niche_breadth = meta$niche_breadth
  )
}


#' Build O(1) occupancy lookups for the carrying-capacity checks.
#'
#' `birth()` recounts occupants with `collapse::fsubset()` on every attempt.
#' Precomputing the counts is equivalent for a frozen state and orders of
#' magnitude faster; `verify_probe()` cross-checks the two against each other.
make_occupancy <- function(agents) {
  ag <- as.data.table(agents)
  inter <- ag[, .(inter_n = .N), by = .(x_loc, y_loc)]
  intra <- ag[, .(intra_n = .N), by = .(x_loc, y_loc, species_id)]
  setkey(inter, x_loc, y_loc)
  setkey(intra, x_loc, y_loc, species_id)
  list(inter = inter, intra = intra)
}


#' Report how much room the capacity ceilings actually leave.
#'
#' If almost every occupied cell sits far below `k_inter`/`k_intra`, the
#' capacity stage is not what separates the modes, and the within-step newborn
#' feedback ignored by the probe cannot matter either.
crowding_headroom <- function(state, k_inter = 50, k_intra = 50) {
  occ <- make_occupancy(state$agents)
  n_hab <- length(state$cells$habitat)
  list(
    occupied_cells = nrow(occ$inter),
    habitat_cells = n_hab,
    mean_inter = mean(occ$inter$inter_n),
    max_inter = max(occ$inter$inter_n),
    frac_cells_at_k_inter = mean(occ$inter$inter_n >= k_inter),
    max_intra = max(occ$intra$intra_n),
    frac_at_k_intra = mean(occ$intra$intra_n >= k_intra)
  )
}


# Core probe ------------------------------------------------------------------

#' Replay birth attempts against a frozen model state.
#'
#' Mirrors the cascade in `birth()` exactly, but records which stage each
#' attempt failed at instead of only returning the survivors.
#'
#' @param state A `full_state` (an element of a recorded `.rds`).
#' @param mode "random" or "short_long", independent of how `state` was
#'   generated - this is what makes the cross-design possible.
#' @param n_trials Number of birth attempts to replay.
#' @param disp Proportion of *short* dispersal, used only when
#'   `mode = "short_long"`. Matches `mod_par$dispersal`.
#'
#' @return A `data.table`, one row per attempt, with the failure stage and the
#'   diagnostic quantities behind it.
establishment_trials <- function(
  state,
  mode = c("random", "short_long"),
  n_trials = 20000L,
  species_par = NULL,
  birth_rate = 0.85,
  k_inter = 50,
  k_intra = 50,
  disp = 1,
  d_dis = NULL,
  kernel = "exponential",
  seed = NULL
) {

  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)

  meta <- state$meta
  grid <- state$grid
  agents <- as.data.table(state$agents)
  habitat_cells <- state$cells$habitat

  if (is.null(species_par)) species_par <- species_par_from_meta(meta)
  if (is.null(d_dis)) d_dis <- meta$dispersal_dist

  grid_size <- nrow(grid)
  nc <- ncol(grid)
  # Values in raster (row-major) cell order, so a destination can be looked up
  # arithmetically instead of through raster's `[` method.
  vals <- raster::getValues(grid)
  nb <- meta$niche_breadth

  occ <- make_occupancy(agents)
  u_by_species <- species_par$n_value  # species_id is the index

  parents <- agents[sample.int(nrow(agents), n_trials, replace = TRUE)]

  out <- data.table(
    trial = seq_len(n_trials),
    mode = mode,
    species_id = parents$species_id,
    u = u_by_species[parents$species_id],
    parent_x = parents$x_loc,
    parent_y = parents$y_loc,
    route = NA_character_,
    dest_x = NA_real_,
    dest_y = NA_real_,
    e = NA_real_,
    survival_prob = NA_real_,
    same_cell = NA,
    stage = NA_character_
  )

  for (i in seq_len(n_trials)) {

    rand <- runif(3, 0, 1)

    if (rand[1] >= birth_rate) {
      set(out, i, "stage", "no_birth")
      next
    }

    cur_loc <- c(parents$x_loc[i], parents$y_loc[i])

    if (mode == "random") {
      route <- "random"
      new_loc <- random_disperse(
        grid = grid,
        habitat_cells = habitat_cells,
        force_habitat = TRUE
      )
    } else {
      if (rand[2] < disp) {
        route <- "short"
        new_loc <- toroidal_disperse(
          cur_loc = cur_loc,
          d_sd = d_dis,
          d_mean = d_dis,
          grid_size = grid_size,
          kernel = kernel
        )
      } else {
        route <- "long"
        new_loc <- random_disperse(grid = grid, force_habitat = FALSE)
      }
    }

    set(out, i, "route", route)
    set(out, i, "dest_x", as.numeric(new_loc[1]))
    set(out, i, "dest_y", as.numeric(new_loc[2]))
    set(out, i, "same_cell", new_loc[1] == cur_loc[1] && new_loc[2] == cur_loc[2])

    e <- vals[(new_loc[1] - 1) * nc + new_loc[2]]

    if (is.na(e)) {
      set(out, i, "stage", "matrix")
      next
    }
    set(out, i, "e", e)

    u <- out$u[i]
    survival_prob <- exp((-(e - u)^2) / (2 * nb^2))
    set(out, i, "survival_prob", survival_prob)

    if (survival_prob <= rand[3]) {
      set(out, i, "stage", "niche")
      next
    }

    inter_n <- occ$inter[.(new_loc[1], new_loc[2]), inter_n]
    intra_n <- occ$intra[.(new_loc[1], new_loc[2], out$species_id[i]), intra_n]
    inter_n <- if (is.na(inter_n)) 0L else inter_n
    intra_n <- if (is.na(intra_n)) 0L else intra_n

    if (inter_n >= k_inter || intra_n >= k_intra) {
      set(out, i, "stage", "capacity")
      next
    }

    set(out, i, "stage", "established")
  }

  out[, stage := factor(
    stage,
    levels = c("no_birth", "matrix", "niche", "capacity", "established")
  )]

  attr(out, "state_label") <- state$step_label
  attr(out, "state_mode") <- meta$dispersal_type
  attr(out, "sim_id") <- meta$sim_id
  out
}


#' Collapse trials into a stage-by-stage conditional decomposition.
summarise_establishment <- function(trials, death_rate = 0.25) {

  n <- nrow(trials)
  n_birth <- trials[stage != "no_birth", .N]
  n_landed <- trials[!stage %in% c("no_birth", "matrix"), .N]
  n_niche_ok <- trials[stage %in% c("capacity", "established"), .N]
  n_est <- trials[stage == "established", .N]

  data.table(
    sim_id = attr(trials, "sim_id"),
    state_mode = attr(trials, "state_mode"),
    probe_mode = trials$mode[1],
    n_trials = n,
    p_birth = n_birth / n,
    p_habitat_given_birth = if (n_birth > 0) n_landed / n_birth else NA_real_,
    p_niche_given_habitat = if (n_landed > 0) n_niche_ok / n_landed else NA_real_,
    p_capacity_given_niche = if (n_niche_ok > 0) n_est / n_niche_ok else NA_real_,
    p_establish = n_est / n,
    death_rate = death_rate,
    net_per_capita = n_est / n - death_rate,
    mean_surv_prob_landed = trials[!is.na(survival_prob), mean(survival_prob)],
    p_same_cell_given_birth = trials[stage != "no_birth", mean(same_cell, na.rm = TRUE)]
  )
}


#' Split short dispersal by whether the offspring landed in the parent's cell.
#'
#' This is the decisive diagnostic. An offspring that lands back in its parent's
#' cell inherits a cell the parent has *already* passed the niche filter on, so
#' its establishment probability is near-certain. `scripts/disp_probs.R` measured
#' that return probability at 24.4% for `mean_disp = 2`. Removing those returns
#' gives the counterfactual "short dispersal with no home-cell advantage", which
#' is the like-for-like comparison against random dispersal.
decompose_by_parent_cell <- function(trials) {
  b <- trials[stage != "no_birth"]
  if (all(is.na(b$same_cell))) stop("No dispersal events recorded.")
  b[, .(
    n = .N,
    share = .N / nrow(b),
    mean_survival_prob = mean(survival_prob, na.rm = TRUE),
    p_established = mean(stage == "established")
  ), by = .(same_cell)][order(same_cell)]
}


#' Counterfactual: how well matched are established agents to their own cells?
#'
#' This is the ceiling that short dispersal approaches as the kernel tightens,
#' and the baseline that random dispersal throws away every generation.
parent_cell_match <- function(state, species_par = NULL) {
  meta <- state$meta
  if (is.null(species_par)) species_par <- species_par_from_meta(meta)
  vals <- raster::getValues(state$grid)
  nc <- ncol(state$grid)
  ag <- as.data.table(state$agents)
  u <- species_par$n_value[ag$species_id]
  e <- vals[(ag$x_loc - 1) * nc + ag$y_loc]
  sp <- exp((-(e - u)^2) / (2 * meta$niche_breadth^2))
  list(
    mean_survival_at_own_cell = mean(sp, na.rm = TRUE),
    median_survival_at_own_cell = median(sp, na.rm = TRUE),
    n_species_present = uniqueN(ag$species_id),
    u_range_present = range(u, na.rm = TRUE)
  )
}


# Verification ----------------------------------------------------------------

#' Cross-check the probe's fast paths against the model's own code.
#'
#' Confirms (1) the arithmetic cell lookup agrees with raster's `[i, j]`, and
#' (2) the precomputed occupancy counts agree with the `collapse::fsubset()`
#' counting that `birth()` actually performs.
verify_probe <- function(state, n_checks = 300L, seed = 1L) {

  set.seed(seed)
  grid <- state$grid
  vals <- raster::getValues(grid)
  nc <- ncol(grid)
  agents <- as.data.table(state$agents)
  occ <- make_occupancy(agents)

  cells <- sample(state$cells$habitat, n_checks, replace = TRUE)
  rc <- cbind(ceiling(cells / nc), (cells - 1) %% nc + 1)

  val_direct <- vapply(seq_len(n_checks), function(k) grid[rc[k, 1], rc[k, 2]], numeric(1))
  val_lookup <- vals[(rc[, 1] - 1) * nc + rc[, 2]]

  spp <- sample(unique(agents$species_id), n_checks, replace = TRUE)
  inter_fsub <- vapply(seq_len(n_checks), function(k) {
    collapse::fnrow(collapse::fsubset(agents, x_loc == rc[k, 1] & y_loc == rc[k, 2]))
  }, numeric(1))
  intra_fsub <- vapply(seq_len(n_checks), function(k) {
    collapse::fnrow(collapse::fsubset(
      agents,
      x_loc == rc[k, 1] & y_loc == rc[k, 2] & species_id == spp[k]
    ))
  }, numeric(1))

  inter_look <- occ$inter[data.table(x_loc = rc[, 1], y_loc = rc[, 2]), inter_n]
  intra_look <- occ$intra[
    data.table(x_loc = rc[, 1], y_loc = rc[, 2], species_id = spp), intra_n
  ]
  inter_look[is.na(inter_look)] <- 0
  intra_look[is.na(intra_look)] <- 0

  list(
    cell_lookup_ok = isTRUE(all.equal(val_direct, val_lookup)),
    inter_count_ok = identical(as.numeric(inter_fsub), as.numeric(inter_look)),
    intra_count_ok = identical(as.numeric(intra_fsub), as.numeric(intra_look))
  )
}


# Analysis --------------------------------------------------------------------

if (sys.nframe() == 0L) {

  # Matched pair: identical ac / hab / frag / edge / nb, differing only in the
  # dispersal mode the simulation was actually run with.
  short_file <- here("output/model_states/0250_ac0_hab0.15_frag0.2_edge1_nb0.1_disp2_r001.rds")
  random_file <- here("output/model_states/0358_ac0_hab0.15_frag0.2_edge1_nb0.1_dispR_r001.rds")

  short_state <- readRDS(short_file)$final
  random_state <- readRDS(random_file)$final

  cat("\n== Probe verification against the model's own code ==\n")
  print(unlist(verify_probe(short_state)))

  cat("\n== Crowding headroom (is the capacity stage binding at all?) ==\n")
  print(unlist(crowding_headroom(short_state)))
  print(unlist(crowding_headroom(random_state)))

  cat("\n== Niche match of agents to their own cells ==\n")
  cat("short-dispersal state:\n")
  print(unlist(parent_cell_match(short_state)))
  cat("random-dispersal state:\n")
  print(unlist(parent_cell_match(random_state)))

  # Cross both probes over both runs at every recorded step. The mode effect is
  # the contrast within a row-pair; the composition effect is the contrast
  # between the two runs for a fixed probe.
  cat("\n== Establishment decomposition across recorded steps ==\n")

  runs <- list(short = short_file, random = random_file)
  results <- rbindlist(lapply(names(runs), function(run) {
    states <- readRDS(runs[[run]])
    rbindlist(lapply(names(states), function(st) {
      rbindlist(lapply(c("short_long", "random"), function(probe) {
        tr <- establishment_trials(states[[st]], mode = probe, n_trials = 20000L, seed = 42)
        cbind(run = run, step = st, summarise_establishment(tr))
      }))
    }))
  }))

  print(results[, .(
    run, step, probe_mode,
    p_hab = round(p_habitat_given_birth, 3),
    p_niche = round(p_niche_given_habitat, 3),
    p_est = round(p_establish, 3),
    net = round(net_per_capita, 3),
    same_cell = round(p_same_cell_given_birth, 3)
  )])

  # The decisive contrast: strip short dispersal's home-cell advantage.
  cat("\n== Short dispersal split by parent-cell return ==\n")
  for (st in c("start", "final")) {
    tr_s <- establishment_trials(
      readRDS(short_file)[[st]], mode = "short_long", n_trials = 30000L, seed = 11
    )
    tr_r <- establishment_trials(
      readRDS(short_file)[[st]], mode = "random", n_trials = 30000L, seed = 11
    )
    cat("\n-- ", st, " state --\n", sep = "")
    print(decompose_by_parent_cell(tr_s), digits = 3)
    away <- tr_s[stage != "no_birth" & same_cell == FALSE]
    cat(sprintf(
      "short, parent-cell returns removed : p_establish = %.3f\nrandom dispersal                   : p_establish = %.3f\n",
      mean(away$stage == "established"),
      mean(tr_r[stage != "no_birth"]$stage == "established")
    ))
  }

  fwrite(results, here("output/establishment_probs.csv"))
  cat("\nWritten to output/establishment_probs.csv\n")
}
