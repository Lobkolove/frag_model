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
#'   `"sp"`) and spatial coordinates (`x_loc`, `y_loc`) for each sampled
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
                 spatvec = NULL,
                 method = c("euclidean", "toroidal", "area"),
                 cutoff = TRUE) {
  
  method <- match.arg(method)
  
  # Drop empty species
  model_sample <- model_sample %>%
    filter(rowSums(across(starts_with("sp"))) > 0)
  
  # Extract species data as presence–absence, drop empty species
  pa_table <- model_sample %>%
    dplyr::select(starts_with("sp")) %>%
    dplyr::select(where(~ sum(.x) > 0)) %>%
    dplyr::mutate(across(everything(), ~ (.x > 0) * 1))
  
  n <- nrow(pa_table)
  
  # Extract coordinates
  coords <- model_sample %>%
    dplyr::select(x_loc, y_loc)

  # Storage for output
  out_list <- vector("list", n)
  
  # Pairwise distances (for euclidean/toroidal)
  if (method %in% c("euclidean", "toroidal")) {
    if (method == "euclidean") {
      pair_dist <- stats::dist(coords) # euclidean distances
    } else if (method == "toroidal") {
      pair_dist <- toroidal_dist(coords, model_sample$grid_size[1]) # toroidal distances 
    }
    pair_dist <- as.matrix(pair_dist)
    
    
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
      S <- rowSums(1 - rich)
      
      spat_ext <- dist_to_site[order(dist_to_site)]
      samp_eff <- seq_len(n)
      
      out_list[[i]] <- dplyr::tibble(
        spat_ext = spat_ext,
        samp_eff = samp_eff,
        S = S
      )
    }
  } else if (method == "area") {
    # Compute full pairwise dist for ordering (euclidean assumed for area ordering)
    pair_dist <- as.matrix(stats::dist(coords))
        
    for (i in 1:n) {
      dist_to_site <- pair_dist[i, ]
      
      # Randomize order to avoid bias for tied distances
      new_order <- sample(seq_len(n))
      new_order <- new_order[order(dist_to_site[new_order])]
      
      # Move focal site to the front
      new_order <- c(i, new_order[new_order != i])
      
      coords_ordered <- as.matrix(coords[new_order, ])
      
      # First point: effort 0
      spat_ext <- numeric(n)
      spat_ext[1] <- 0
      
      # Second point: distance to first (non-zero area starts at 3rd)
      spat_ext[2] <- dist_to_site[order(dist_to_site)][1]  # or 0, but match mobr
      
      # Convex hull areas for 3 to n
      for (j in 3:n) {
        hpts <- grDevices::chull(coords_ordered[1:j, 1], coords_ordered[1:j, 2])
        hpts <- c(hpts, hpts[1])
        chull_coords <- as.matrix(coords_ordered[hpts, ])
        chull_poly <- sf::st_polygon(list(chull_coords))
        chull_area <- sf::st_area(chull_poly)
        spat_ext[j] <- as.numeric(chull_area)
      }

      comm_ordered <- pa_table[new_order, ]
      comm_bool <- (comm_ordered == 0) * 1
      rich <- apply(comm_bool, 2, cumprod)
      S <- rowSums(1 - rich)
      samp_eff <- seq_len(n)

      out_list[[i]] <- dplyr::tibble(
        spat_ext = spat_ext,
        samp_eff = samp_eff,
        S = S
      )
    }
  }
  
  # Compute richness curves (common to all)
  scr_mat <- matrix(0, n, n)
  for (i in 1:n) {
    new_order <- sample(seq_len(n))
    new_order <- new_order[order(pair_dist[i, new_order])]
    new_order <- c(i, new_order[new_order != i])
    
    comm_ordered <- pa_table[new_order, ]
    comm_bool <- (comm_ordered == 0) * 1
    rich <- apply(comm_bool, 2, cumprod)
    
    scr_mat[i, ] <- ncol(pa_table) - rowSums(rich)
  }
  
  # Combine output into a single data frame
  out_dat <- dplyr::bind_rows(out_list)
  
  # Cutoff logic
  if (cutoff & method == "euclidean") {
    d_cut <- model_sample$grid_size[1] / 2
    out_dat <- out_dat %>%
      dplyr::filter(spat_ext <= d_cut)
  } else if (cutoff & method == "area") {
    gs <- model_sample$grid_size[1]
    a_cut <- (gs ^ 2) / 4  # 1/4 grid surface
    out_dat <- out_dat %>%
      dplyr::filter(spat_ext <= a_cut)
  }

  # Fit model - GAM with monotonously increasing constraint
  scam1 <- scam::scam(S ~ s(spat_ext, bs = "mpi") + s(samp_eff),
                      data = out_dat, family = "poisson")
  
  # Effort grid for interpolation
  if (is.null(spatvec)) {
    spatvec <- seq(min(out_dat$spat_ext, na.rm = TRUE),
                   max(out_dat$spat_ext, na.rm = TRUE),
                   length = 200)
  }
  
  # Predict sampling effort for each spatial extent
  lm_ss <- lm(samp_eff ~ spat_ext, data = out_dat)


  samp_eff_pred <- predict(lm_ss, newdata = data.frame(spat_ext = spatvec))

  out_pred <- data.frame(spat_ext = spatvec, samp_eff = samp_eff_pred)
  
  pred <- predict(scam1, out_pred, se = TRUE, type = "response")

  out_pred$S <- pred$fit
  out_pred$S_low <- pred$fit - 2*pred$se.fit
  out_pred$S_high <- pred$fit + 2*pred$se.fit

  out_pred <- out_pred %>%
    group_by(spat_ext) %>%
    summarise(S = mean(S), 
              S_low = mean(S_low), 
              S_high = mean(S_high))
  
  out <- list(data   = out_dat,    
              smooth = out_pred)   
  
  class(out) <- "sSBR"
  attr(out, "method") <- method
  return(out)
}

plot.sSBR <- function(sSBR_object,
                      col = "midnightblue",
                      all_lines = TRUE,
                      ...) {
  
  dat <- sSBR_object$data
  sm  <- sSBR_object$smooth
  meth <- attr(sSBR_object, "method")
  
  # Dynamic x-label based on method
  xlab <- switch(meth,
                 "euclidean" = "Euclidean distance",
                 "toroidal"  = "Toroidal distance",
                 "area"      = "Cumulative convex hull area")
  
  # Transparency scales with number of curves (only if plotting them)
  if (all_lines) {
    a <- 1 / sqrt(nrow(dat) / 300)
    alpha <- pmax(0.05, pmin(0.75, a))
  }
  
  # Base plot setup
  graphics::plot(NA,
                 xlim = range(dat$effort, na.rm = TRUE),
                 ylim = range(dat$S, na.rm = TRUE),
                 xlab = xlab,
                 ylab = "Cumulative species richness",
                 las = 1)
  
  # Individual rarefaction curves (optional, semi-transparent)
  if (all_lines) {
    for (i in unique(dat$id)) {
      tmp <- dat[dat$id == i, ]
      graphics::lines(tmp$effort,
                      tmp$S,
                      col = colorspace::adjust_transparency(col, alpha = alpha))
    }
  }
  
  # 95% CI ribbon for smooth
  graphics::polygon(c(sm$effort, rev(sm$effort)),
                    c(sm$S_low, rev(sm$S_high)),
                    col = colorspace::adjust_transparency(col, alpha = .5),
                    border = NA)
  
  # Central smooth line (always plotted, bold)
  graphics::lines(sm$effort,
                  sm$S,
                  col = colorspace::darken(col, amount = .4),
                  lwd = 3)
  
  invisible(sSBR_object)
}

