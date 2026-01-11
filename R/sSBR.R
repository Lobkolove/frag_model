#' Spatially constrained rarefaction from model sample output
#'
#' Computes a spatially explicit, sample-based species rarefaction curve
#' (sSBR sensu McGlinn et al. 2019), also referred to as a spatially
#' constrained rarefaction (SCR; Chiarucci et al. 2009), from
#' sample-level output of the agent-based model.
#'
#' The function is tailored to model outputs representing spatially
#' explicit community samples with known locations and species
#' abundances, and is not applicable to aggregated or time-step–level
#' summary outputs.
#'
#' @param model_sample A data frame corresponding to *sample-level*
#'   model output, containing species abundance columns (prefixed with
#'   `"sp"`) and spatial coordinates (`loc_x`, `loc_y`) for each sampled
#'   location.
#' @param distvec Optional numeric vector of spatial distances at which
#'   the rarefaction curve should be interpolated. If `NULL`, a regular
#'   sequence spanning the observed distance range is used.
#'
#' @details
#' Species abundances are converted to presence–absence data and species
#' with no occurrences across all samples are removed prior to analysis.
#' For each focal sample, neighboring samples are ordered by increasing
#' spatial distance, with random permutation applied to break ties and
#' avoid ordering bias.
#'
#' The curve is expressed as a function of accumulated spatial distance
#' rather than accumulated sample count. A monotonically increasing
#' generalized additive model (GAM) is used to interpolate cumulative
#' species richness as a function of distance.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{data}{A data frame containing cumulative spatial
#'       distances and corresponding cumulative species richness values
#'       for all focal samples.}
#'     \item{smooth}{A data frame containing the interpolated SCR
#'       curve and approximate confidence intervals derived from a
#'       monotonic GAM fit.}
#'   }
#'
#' @references
#' Chiarucci, A. et al. (2009). Spatially constrained rarefaction:
#' incorporating the autocorrelated structure of biological communities
#' into sample-based rarefaction. *Community Ecology*, 10, 209–214.
#'
#' McGlinn, D. J. et al. (2019). Measurement of Biodiversity (MoB):
#' A method to separate the scale-dependent effects of species abundance
#' distribution, density, and aggregation on diversity change.
#' *Methods in Ecology and Evolution*, 10, 258–269.
#'
#' @export
#' 
#' @noRd

sSBR <- function(model_sample, 
                 distvec = NULL) {
  
  # Extract species data as presence–absence, drop empty species
  pa_table <- model_sample %>%
    select(starts_with("sp")) %>%
    select(where(~ sum(.x) > 0)) %>%
    mutate(across(everything(), ~ as.integer(.x > 0)))
  
  n <- nrow(pa_table)
  
  # Extract coordinates
  coords <- model_sample %>%
    filter(rowSums(across(starts_with("sp"))) > 0) %>%
    select(loc_x, loc_y)
  
  # Pairwise distances
  pair_dist <- as.matrix(stats::dist(coords))
  
  # Storage
  scr_mat  <- matrix(0, n, n)
  dist_mat <- matrix(0, n, n)
  
  for (i in 1:n) {
    
    dist_to_site <- pair_dist[i, ]
    
    # Randomize order to avoid bias for tied distances
    new_order <- sample(seq_len(n))
    new_order <- new_order[order(dist_to_site[new_order])]
    
    # Move focal site to the front
    new_order <- c(i, new_order[new_order != i])
    
    comm_ordered <- pa_table[new_order, ]
    
    # 1 for absence, 0 for presence
    comm_bool <- (comm_ordered == 0) * 1
    rich <- apply(comm_bool, 2, cumprod)

    
    scr_mat[i, ]  <- ncol(pa_table) - rowSums(rich)
    dist_mat[i, ] <- dist_to_site[order(dist_to_site)]
  }
  
  # Long-format data: distance–richness pairs
  out_dat <- data.frame(id       = rep(1:n, times = n),
                        distance = as.vector(dist_mat),
                        S        = as.vector(scr_mat))
  
  out_dat <- out_dat[order(out_dat$id, out_dat$distance), ]
  
  
  # Fit model - GAM with monotonously increasing constraint
  scam1 <- scam::scam(S ~ s(distance, bs = "mpi"),
                      data = out_dat, family = "poisson")
  
  # Distance grid for interpolation
  if (is.null(distvec)) {
    distvec <- seq(min(out_dat$distance),
                   max(out_dat$distance),
                   length = 200)
  }
  
  out_pred <- data.frame(distance = distvec, 
                         S = NA)
  
  pred <- predict(scam1, out_pred, se = TRUE, type = "response")
  
  out_pred$S <- pred$fit
  out_pred$S_low <- pred$fit - 2*pred$se.fit
  out_pred$S_high <- pred$fit + 2*pred$se.fit
  # These confidence intervals ignore the dependence between the points,
  # so they are likely inappropiate.
  
  out <- list(data   = out_dat,    # id, distance, S
              smooth = out_pred)    # distance, S, CI
  
  class(out) <- "sSBR"
  return(out)
  
}

plot.sSBR <- function(sSBR_object,
                      col = "blue",
                      ...) {
  
  dat <- sSBR_object$data
  sm  <- sSBR_object$smooth
  
  # Base plot
  graphics::plot(NA,
                 xlim = range(dat$distance, na.rm = TRUE),
                 ylim = range(dat$S, na.rm = TRUE),
                 xlab = "Spatial distance",
                 ylab = "Cumulative species richness",
                 las = 1,
                 main = "Spatially constrained rarefaction")
  
  # Individual curves (one line per sample)
  for (i in unique(dat$id)) {
    tmp <- dat[dat$id == i, ]
    graphics::lines(tmp$distance,
                    tmp$S,
                    col = adjustcolor(col, alpha.f = .5))
  }
  
  # Confidence ribbon
  graphics::polygon(c(sm$distance, rev(sm$distance)),
                    c(sm$S_low, rev(sm$S_high)),
                    col = adjustcolor(col, alpha.f = .35),
                    border = NA)

  
  # Prediction line
  graphics::lines(sm$distance,
                  sm$S,
                  col = col,
                  lwd = 3)
  
  invisible(sSBR_object)
}



