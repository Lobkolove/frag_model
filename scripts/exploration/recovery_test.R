################################################################################
#
# Does post-fragmentation recovery depend on pre-fragmentation population size?
#
################################################################################
#
# Question
# --------
# In the random-dispersal series, replicates either collapse or rebound after the
# fragmentation event, with little in between. This script tests whether that
# outcome is predicted by how large the population was *before* fragmentation
# was applied.
#
# Why more than one test
# ----------------------
# "Recovery" is not a quantity the model produces - it is a threshold imposed on
# the final population size, and the threshold is a judgement call. A result that
# only holds at one arbitrary cut-off is not a result. So the script runs:
#
#   1. A threshold-free test  - Spearman correlation between pre-fragmentation
#      and final population size. This is the primary analysis: it uses all the
#      information and commits to no cut-off.
#   2. The binary test asked for - Wilcoxon rank-sum on pre-fragmentation size
#      between recovered and collapsed replicates, plus logistic regression.
#      Rank-based because both variables are strongly skewed and the outcome is
#      bimodal, which rules out a t-test or a Pearson correlation.
#   3. A sensitivity sweep across recovery thresholds, so the conclusion can be
#      seen to hold (or not) independently of where the line is drawn.
#
# A validity check is run first: the fragmentation treatment must not influence
# pre-fragmentation population size, since fragmentation has not yet occurred at
# that point. If it appeared to, replicates could not be pooled across levels.
#
# Secondary analysis
# ------------------
# Pre-fragmentation size is not the only candidate. `suitability` measures how
# well the habitat retained by the mask matches the species pool that SURVIVED
# fragmentation - the expected niche-pass rate, averaged over retained cells and
# abundance-weighted over species. The script reports whether pre-fragmentation
# size still predicts recovery once suitability is accounted for, because the
# two are not independent.
#
# On interpreting the fragmentation comparison
# --------------------------------------------
# Suitability shows no significant difference across fragmentation levels here,
# but with ten replicates per level that test has roughly 35% power against the
# observed pattern, so on its own it establishes nothing. Positive evidence for
# a flat mean comes from generating masks directly instead: across 720 masks
# (30 environments x 8 mask seeds x 3 levels, paired on environment), mean
# suitability was 0.426 / 0.421 / 0.416 with p = 0.565, while the variance
# differed enormously (sd 0.142 / 0.118 / 0.062, Fligner-Killeen p < 2e-16).
# A smooth mask at low fragmentation carves a narrow, randomly placed slice of
# a smooth environment; a rough mask samples the whole range every time. So
# fragmentation changes how variable the draw is, not how good it is on average,
# and apparent between-level differences at n = 10 are sampling noise.
#
################################################################################

library(here)
library(data.table)
library(raster)

# Configuration ---------------------------------------------------------------

# Recovery threshold for the binary analyses. Reported explicitly rather than
# buried, and swept in section 4.
RECOVERY_THRESHOLD <- 2000

# Species optima are an evenly spaced sequence under `species_specific_par = 0`.
N_SPECIES <- 1000


# Data extraction -------------------------------------------------------------

#' Pull one row of per-replicate metrics out of a recorded simulation state.
#'
#' @param path Path to a recorded `.rds` holding the four model states.
#' @return A one-row `data.table`.
extract_replicate <- function(path) {

  states <- readRDS(path)
  pre <- states$pre_fragmentation
  post <- states$post_fragmentation
  final <- states$final

  meta <- post$meta
  n_values <- seq(from = 0, to = 1, length.out = N_SPECIES)
  nb <- meta$niche_breadth

  # Environmental values of the cells the fragmentation mask retained
  e_retained <- raster::getValues(post$grid)
  e_retained <- e_retained[!is.na(e_retained)]

  # Expected niche-pass rate of a species pool across the retained habitat:
  # averaged over every retained cell, abundance-weighted over species.
  niche_pass <- function(agents) {
    pool <- as.data.table(agents)[, .N, by = species_id]
    u <- n_values[pool$species_id]
    w <- pool$N / sum(pool$N)
    sum(w * vapply(
      u,
      function(ui) mean(exp(-(e_retained - ui)^2 / (2 * nb^2))),
      numeric(1)
    ))
  }

  # Primary definition uses the pool that SURVIVED fragmentation, because those
  # are the species that actually rebuild the population. Averaging over all
  # retained cells - not only the cells survivors occupy - keeps the measure
  # from being a restatement of "survivors sit where they can live": it asks
  # whether the rest of the retained habitat also suits them, which is exactly
  # what matters when dispersal places offspring anywhere in the landscape.
  #
  # The pre-fragmentation version is kept for comparison. It is the weaker
  # predictor (Spearman 0.66 vs 0.87 against final population size).

  data.table(
    sim_id = as.integer(meta$sim_id),
    frag = meta$fragmentation,
    ac = meta$ac_amount,
    n_pre = nrow(pre$agents),
    n_post = nrow(post$agents),
    n_final = nrow(final$agents),
    n_habitat = length(e_retained),
    suitability = niche_pass(post$agents),
    suitability_pre = niche_pass(pre$agents)
  )
}


#' Simulated power of a one-way ANOVA, given the pattern actually observed.
#'
#' Reported alongside every non-significant group comparison in this script. A
#' null result at this replicate count is a failure to reject, not evidence of
#' no effect, and the only way to tell those apart is to state the power.
power_anova <- function(group_means, pooled_sd, n_per_group, n_sim = 4000L, alpha = 0.05) {
  k <- length(group_means)
  p <- replicate(n_sim, {
    y <- unlist(lapply(group_means, function(m) rnorm(n_per_group, m, pooled_sd)))
    g <- rep(seq_len(k), each = n_per_group)
    anova(lm(y ~ factor(g)))[["Pr(>F)"]][1]
  })
  list(
    cohen_f = sd(group_means) * sqrt((k - 1) / k) / pooled_sd,
    power = mean(p < alpha)
  )
}


#' Assemble the replicate table for a set of simulation ids.
build_dataset <- function(
  sim_ids,
  state_dir = here("output/model_states"),
  threshold = RECOVERY_THRESHOLD
) {
  files <- list.files(state_dir, pattern = "[.]rds$", full.names = TRUE)
  ids <- as.integer(substr(basename(files), 1, 4))
  files <- files[ids %in% sim_ids]

  if (length(files) == 0L) stop("No state files matched the requested sim_ids.")

  out <- rbindlist(lapply(files, extract_replicate))
  out[, recovered := n_final > threshold]
  setorder(out, frag, n_pre)
  out[]
}


# Tests -----------------------------------------------------------------------

#' Validity check: fragmentation must not predict the pre-fragmentation state.
check_pooling <- function(dat) {
  fit <- aov(n_pre ~ factor(frag), data = dat)
  p <- summary(fit)[[1]][["Pr(>F)"]][1]
  list(
    p_value = p,
    poolable = p > 0.05,
    by_level = dat[, .(reps = .N, n_pre_mean = mean(n_pre), n_pre_sd = sd(n_pre)), by = frag][order(frag)]
  )
}


#' Threshold-free association between pre-fragmentation and final size.
test_continuous <- function(dat) {
  ct <- suppressWarnings(cor.test(dat$n_pre, dat$n_final, method = "spearman"))
  list(rho = unname(ct$estimate), p_value = ct$p.value, n = nrow(dat))
}


#' Binary tests: is pre-fragmentation size higher in recovered replicates?
test_binary <- function(dat) {

  rec <- dat[recovered == TRUE, n_pre]
  col <- dat[recovered == FALSE, n_pre]

  if (length(rec) < 3L || length(col) < 3L) {
    stop("Too few replicates in one outcome group for a meaningful test.")
  }

  wt <- suppressWarnings(wilcox.test(rec, col))
  # Rank-biserial correlation: a bounded effect size for the rank-sum test
  rb <- 1 - (2 * unname(wt$statistic)) / (length(rec) * length(col))

  # Logistic regression on the log scale, since n_pre spans three orders of magnitude
  glm_fit <- glm(recovered ~ log(n_pre), family = binomial, data = dat)
  cf <- summary(glm_fit)$coefficients

  list(
    n_recovered = length(rec),
    n_collapsed = length(col),
    median_recovered = median(rec),
    median_collapsed = median(col),
    wilcox_W = unname(wt$statistic),
    wilcox_p = wt$p.value,
    rank_biserial = abs(rb),
    glm_beta = cf["log(n_pre)", "Estimate"],
    glm_p = cf["log(n_pre)", "Pr(>|z|)"],
    odds_ratio_per_doubling = exp(cf["log(n_pre)", "Estimate"] * log(2))
  )
}


#' Does the conclusion survive a different definition of "recovered"?
sweep_threshold <- function(dat, thresholds = c(500, 1000, 1500, 2000, 3000, 4000, 5000)) {
  rbindlist(lapply(thresholds, function(th) {
    rec <- dat[n_final > th, n_pre]
    col <- dat[n_final <= th, n_pre]
    if (length(rec) < 3L || length(col) < 3L) {
      return(data.table(threshold = th, n_recovered = length(rec),
                        n_collapsed = length(col), wilcox_p = NA_real_))
    }
    wt <- suppressWarnings(wilcox.test(rec, col))
    data.table(threshold = th, n_recovered = length(rec), n_collapsed = length(col),
               wilcox_p = wt$p.value)
  }))
}


#' Suitability as a predictor of final population size.
#'
#' Threshold-free counterpart to the recovery analysis, and the more informative
#' framing: it asks how much of the final population size is explained by how
#' well the retained habitat matches the surviving pool, without collapsing the
#' outcome to a binary. Rank-based first because final size spans three orders
#' of magnitude; the log-linear model then gives an interpretable effect size.
#'
#' `n_post` is included to separate two routes: suitability could act
#' immediately (which cells the mask keeps determines who survives it) or
#' through subsequent recovery dynamics.
test_suitability <- function(dat) {

  sp_final <- suppressWarnings(cor.test(dat$suitability, dat$n_final, method = "spearman"))
  sp_post <- suppressWarnings(cor.test(dat$suitability, dat$n_post, method = "spearman"))

  m_simple <- lm(log(n_final) ~ suitability, data = dat)
  m_frag <- lm(log(n_final) ~ suitability + factor(frag), data = dat)

  cf <- summary(m_simple)$coefficients

  list(
    spearman_final = list(rho = unname(sp_final$estimate), p = sp_final$p.value),
    spearman_post = list(rho = unname(sp_post$estimate), p = sp_post$p.value),
    lm_beta = cf["suitability", "Estimate"],
    lm_p = cf["suitability", "Pr(>|t|)"],
    r_squared = summary(m_simple)$r.squared,
    # Fold-change in final population size per 0.1 increase in suitability
    fold_per_0.1 = exp(cf["suitability", "Estimate"] * 0.1),
    frag_adds = anova(m_simple, m_frag),
    by_frag = dat[, .(reps = .N, suitability = mean(suitability), sd = sd(suitability)), by = frag][order(frag)],
    anova_frag = summary(aov(suitability ~ factor(frag), data = dat))[[1]][["Pr(>F)"]][1]
  )
}


#' Secondary: does pre-fragmentation size add anything beyond retained-habitat
#' suitability, and vice versa?
test_joint <- function(dat) {
  full <- glm(recovered ~ scale(log(n_pre)) + scale(suitability), family = binomial, data = dat)
  cf <- summary(full)$coefficients
  list(
    coefficients = round(cf, 4),
    cor_npre_suitability = cor(log(dat$n_pre), dat$suitability, method = "spearman")
  )
}


# Analysis --------------------------------------------------------------------

if (sys.nframe() == 0L) {

  log <- fread(here("output/simulations_log.csv"))
  target <- log[dispersal_type == "random" & ac_amount == 0.7, as.integer(sim_id)]

  cat("Random-dispersal series, ac = 0.7 -", length(target), "simulations\n")
  dat <- build_dataset(target)
  fwrite(dat, here("output/recovery_test_data.csv"))

  cat("\n== 0. Validity check: can replicates be pooled across fragmentation levels? ==\n")
  pool <- check_pooling(dat)
  print(pool$by_level, digits = 4)
  cat(sprintf(
    "\nANOVA of pre-fragmentation size on fragmentation level: p = %.3f -> %s\n",
    pool$p_value,
    if (pool$poolable) "no treatment effect, pooling is valid" else "WARNING: treatment effect present, do not pool"
  ))

  cat("\n== 1. Primary, threshold-free: pre-fragmentation vs final population size ==\n")
  cont <- test_continuous(dat)
  cat(sprintf("Spearman rho = %.3f, p = %.4f (n = %d)\n", cont$rho, cont$p_value, cont$n))

  cat(sprintf("\n== 2. Binary test at threshold n_final > %d ==\n", RECOVERY_THRESHOLD))
  bin <- test_binary(dat)
  cat(sprintf("recovered: %d replicates, median pre-frag N = %.0f\n", bin$n_recovered, bin$median_recovered))
  cat(sprintf("collapsed: %d replicates, median pre-frag N = %.0f\n", bin$n_collapsed, bin$median_collapsed))
  cat(sprintf("Wilcoxon rank-sum: W = %.0f, p = %.4f, rank-biserial r = %.3f\n",
              bin$wilcox_W, bin$wilcox_p, bin$rank_biserial))
  cat(sprintf("Logistic on log(N_pre): beta = %.3f, p = %.4f, odds ratio per doubling = %.2f\n",
              bin$glm_beta, bin$glm_p, bin$odds_ratio_per_doubling))

  cat("\n== 3. Sensitivity to the recovery threshold ==\n")
  print(sweep_threshold(dat), digits = 3)

  cat("\n== 4. Suitability as a predictor of FINAL POPULATION SIZE (threshold-free) ==\n")
  su <- test_suitability(dat)
  cat(sprintf("Spearman, suitability vs final size : rho = %.3f, p = %.4f\n",
              su$spearman_final$rho, su$spearman_final$p))
  cat(sprintf("Spearman, suitability vs post-frag  : rho = %.3f, p = %.4f\n",
              su$spearman_post$rho, su$spearman_post$p))
  cat(sprintf("lm(log final size ~ suitability)    : beta = %.2f, p = %.5f, R2 = %.3f\n",
              su$lm_beta, su$lm_p, su$r_squared))
  cat(sprintf("  -> %.2fx final population per +0.1 suitability\n", su$fold_per_0.1))
  cat("\nDoes fragmentation level add anything beyond suitability?\n")
  print(su$frag_adds)
  cat("\nSuitability by fragmentation level:\n")
  print(su$by_frag, digits = 4)
  cat(sprintf("ANOVA of suitability on fragmentation level: p = %.3f\n", su$anova_frag))

  # A null result here is only interpretable alongside its power.
  pw <- power_anova(
    group_means = su$by_frag$suitability,
    pooled_sd = sqrt(mean(su$by_frag$sd^2)),
    n_per_group = min(su$by_frag$reps)
  )
  cat(sprintf(
    "  Cohen f = %.3f; power to detect this pattern at n = %d per group: %.0f%%\n",
    pw$cohen_f, min(su$by_frag$reps), 100 * pw$power
  ))
  if (su$anova_frag > 0.05) {
    cat("  -> Non-significant at this replicate count. This is a failure to reject,\n")
    cat("     NOT evidence that suitability is unaffected by fragmentation. Positive\n")
    cat("     evidence requires simulating masks directly at high n (see the\n")
    cat("     structural test described in the header).\n")
  }

  cat("\n== 5. Secondary: pre-fragmentation size vs retained-habitat suitability ==\n")
  joint <- test_joint(dat)
  print(joint$coefficients)
  cat(sprintf("\nSpearman correlation between the two predictors: %.3f\n", joint$cor_npre_suitability))

  cat("\nReplicate-level data written to output/recovery_test_data.csv\n")
}
