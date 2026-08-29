################################################################################
#
# Population and richness trends across the random-dispersal series
#
################################################################################
#
# Four states are recorded per simulation, at steps 0, 40, 40 and 100. Both
# fragmentation states share step 40 - `fragment()` does not advance the step
# counter - so the fragmentation event appears as a vertical drop at the dashed
# line. This is why the curves are drawn with `geom_path()` rather than
# `geom_line()`: geom_line sorts points by x and would connect the two step-40
# states in whichever order it liked.
#
# Replicate-level paths are drawn faintly behind the group means. This is not
# decoration: final population size in this series is strongly bimodal, with
# replicates either rebounding to several thousand individuals or staying near
# extinction. A group mean sits in the gap between those two outcomes and
# describes almost none of the replicates, so it should not be read alone.
#
################################################################################

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

theme_set(theme_bw(base_size = 14))

# Palette shared with scripts/clean_diversity.R
pal_frag <- c("#4688ad", "#e9b14a", "#ac5384")

FRAG_STEP <- 40


# Data ------------------------------------------------------------------------

#' Per-simulation, per-step abundance and richness for a set of simulations.
#'
#' The `fragmentation` column in the exported data is NA for the `start` and
#' `pre_fragmentation` states, because `clean_run()` only writes the value into
#' the metadata once the fragmentation event has happened. Grouping on that
#' column directly would therefore split every simulation's early steps into an
#' NA level and break the trend lines. The level each simulation was run at is
#' taken from the log instead, and applied to all four of its states.
collect_trends <- function(
  sim_ids,
  sampled_dir = here("output/sampled_data"),
  log_file = here("output/simulations_log.csv")
) {

  files <- list.files(sampled_dir, pattern = "_all[.]csv$", full.names = TRUE)
  ids <- as.integer(substr(basename(files), 1, 4))
  files <- files[ids %in% sim_ids]
  if (length(files) == 0L) stop("No sampled files matched the requested sim_ids.")

  frag_lookup <- fread(log_file)[, .(sim_id = as.integer(sim_id), frag_level = fragmentation)]

  out <- rbindlist(lapply(files, function(p) {
    d <- fread(p)[!is.na(species_id) & n > 0]
    d[, .(N = sum(n), SR = uniqueN(species_id)),
      by = .(sim_id, step, step_label)]
  }))

  out <- merge(out, frag_lookup, by = "sim_id", all.x = TRUE)
  if (anyNA(out$frag_level)) stop("Some simulations are missing a fragmentation level in the log.")
  setnames(out, "frag_level", "fragmentation")

  # Explicit ordering so the two step-40 states are drawn pre- then post-
  out[, step_order := match(
    step_label,
    c("start", "pre_fragmentation", "post_fragmentation", "final")
  )]
  out[, fragmentation := factor(
    fragmentation,
    levels = c(0.2, 0.5, 0.8),
    labels = c("Low", "Medium", "High")
  )]
  setorder(out, sim_id, step_order)
  out[]
}


# Plot ------------------------------------------------------------------------

#' Trend plot for one response variable.
#'
#' @param trends Output of `collect_trends()`.
#' @param var "N" or "SR".
#' @param log_y Log-scale the y axis. Final population size spans 20 to ~13,000
#'   in this series, so on a linear axis the collapsed replicates are pinned to
#'   the baseline and invisible.
#' @param show_replicates Draw individual replicate paths behind the means.
plot_trend <- function(trends, var = c("N", "SR"), log_y = FALSE, show_replicates = TRUE) {

  var <- match.arg(var)
  y_lab <- c(N = "Number of individuals", SR = "Species richness")[[var]]

  dat <- copy(trends)
  dat[, y := get(var)]

  means <- dat[, .(y = mean(y)), by = .(fragmentation, step, step_order, step_label)]
  setorder(means, fragmentation, step_order)

  p <- ggplot(mapping = aes(x = step, y = y, colour = fragmentation)) +
    geom_vline(
      xintercept = FRAG_STEP,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.6
    )

  if (show_replicates) {
    p <- p + geom_path(
      data = dat,
      aes(group = sim_id),
      alpha = 0.18,
      linewidth = 0.4
    )
  }

  p <- p +
    geom_path(data = means, aes(group = fragmentation), linewidth = 1.2) +
    geom_point(data = means, size = 2.6) +
    scale_colour_manual(values = pal_frag) +
    scale_x_continuous(breaks = c(0, FRAG_STEP, 100)) +
    labs(
      x = "Time step",
      y = y_lab,
      colour = "Level of\nfragmentation"
    )

  if (log_y) {
    p <- p + scale_y_log10() + labs(y = paste0(y_lab, " (log scale)"))
  }

  p
}


#' Distribution of population size in the step before fragmentation.
#'
#' Drawn on a log x axis because the values span 327 to 18,262. Fill is by
#' fragmentation level purely to show that the levels are interleaved rather
#' than separated - fragmentation has not been applied at this point, so any
#' apparent grouping here would indicate a problem with the design.
plot_pre_frag_spread <- function(trends, bins = 14) {

  pre <- trends[step_label == "pre_fragmentation"]

  ggplot(pre, aes(x = N)) +
    geom_histogram(
      aes(fill = fragmentation),
      bins = bins,
      colour = "white",
      linewidth = 0.3
    ) +
    geom_vline(
      xintercept = 5000,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.6
    ) +
    geom_vline(
      xintercept = median(pre$N),
      linetype = "solid",
      colour = "grey20",
      linewidth = 0.6
    ) +
    geom_rug(colour = "grey30", alpha = 0.7) +
    annotate(
      "text", x = 5000, y = Inf, label = "  starting N", hjust = 0, vjust = 1.6,
      colour = "grey40", size = 3.4
    ) +
    annotate(
      "text", x = median(pre$N), y = Inf, label = "median  ", hjust = 1, vjust = 1.6,
      colour = "grey20", size = 3.4
    ) +
    scale_x_log10(breaks = c(300, 1000, 3000, 10000)) +
    scale_fill_manual(values = pal_frag) +
    labs(
      x = "Number of individuals before fragmentation (log scale)",
      y = "Replicates",
      fill = "Level of\nfragmentation"
    )
}


# Output ----------------------------------------------------------------------

if (sys.nframe() == 0L) {

  log <- fread(here("output/simulations_log.csv"))
  target <- log[dispersal_type == "random" & ac_amount == 0.7, as.integer(sim_id)]
  cat("Random-dispersal series, ac = 0.7 -", length(target), "simulations\n")

  trends <- collect_trends(target)

  gg_n <- plot_trend(trends, "N")
  gg_sr <- plot_trend(trends, "SR")
  gg_n_log <- plot_trend(trends, "N", log_y = TRUE)
  gg_sr_log <- plot_trend(trends, "SR", log_y = TRUE)

  ggsave(here("pics/random_trends_N.png"), gg_n, width = 7, height = 5, dpi = 300)
  ggsave(here("pics/random_trends_SR.png"), gg_sr, width = 7, height = 5, dpi = 300)
  ggsave(here("pics/random_trends_N_log.png"), gg_n_log, width = 7, height = 5, dpi = 300)
  ggsave(here("pics/random_trends_SR_log.png"), gg_sr_log, width = 7, height = 5, dpi = 300)

  gg_combined <- (gg_n_log / gg_sr_log) +
    plot_layout(guides = "collect") &
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 14), legend.position = "right")

  ggsave(here("pics/random_trends_combined.png"), gg_combined, width = 8, height = 9, dpi = 300)

  gg_spread <- plot_pre_frag_spread(trends)
  ggsave(here("pics/random_pre_frag_spread.png"), gg_spread, width = 7, height = 4.5, dpi = 300)

  pre <- trends[step_label == "pre_fragmentation"]
  cat("\nPre-fragmentation population size, all replicates:\n")
  print(sort(pre$N))
  cat(sprintf(
    "\nmin %d | Q1 %.0f | median %.0f | mean %.0f | Q3 %.0f | max %d | sd %.0f | CV %.2f\n",
    min(pre$N), quantile(pre$N, .25), median(pre$N), mean(pre$N),
    quantile(pre$N, .75), max(pre$N), sd(pre$N), sd(pre$N) / mean(pre$N)
  ))
  cat(sprintf("below the starting 5000: %d of %d replicates\n", sum(pre$N < 5000), nrow(pre)))

  cat("\nGroup means by fragmentation level:\n")
  print(dcast(
    trends[, .(N = round(mean(N)), SR = round(mean(SR))), by = .(fragmentation, step_label)],
    step_label ~ fragmentation,
    value.var = c("N", "SR")
  ))

  cat("\nWritten to pics/random_trends_{N,SR}.png (linear),",
      "\n            pics/random_trends_{N,SR}_log.png (log y),",
      "\n            pics/random_trends_combined.png\n")
}
