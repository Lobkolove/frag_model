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
                 method = c("euclidean", "toroidal", "area"),
                 cutoff = TRUE,
                 spatvec = NULL,
                 effort_ref = NULL,
                 n_focal = NULL) {
  
  method <- match.arg(method)
  
  # Drop empty species rows
  model_sample <- model_sample %>%
    dplyr::filter(rowSums(dplyr::across(starts_with("sp"))) > 0)
  
  # Shared spatial accumulation + richness computation
  out_dat <- spat_acc(
    sampled_data = model_sample,
    method = method,
    n_focal = n_focal,
    compute_richness = TRUE
  )
  
  # Cutoff logic
  if (cutoff & method == "euclidean") {
    d_cut <- model_sample$grid_size[1] / 2
    out_dat <- out_dat %>% dplyr::filter(spat_ext <= d_cut)
  } else if (cutoff & method == "area") {
    gs <- model_sample$grid_size[1]
    a_cut <- (gs ^ 2) / 4  # 1/4 grid surface
    out_dat <- out_dat %>% dplyr::filter(spat_ext <= a_cut)
  }
  
  # Fit model - GAM with monotonously increasing constraint
  scam1 <- scam::scam(S ~ s(spat_ext, bs = "mpi") + s(samp_eff),
                      data = out_dat, family = "poisson")
  
  # Prepare prediction grid
  if (!is.null(effort_ref)) {
    # Use provided standardized effort reference
    out_pred <- effort_ref
  } else {
    # Fallback: scenario-specific median effort per spat_ext
    if (is.null(spatvec)) {
      spatvec <- seq(min(out_dat$spat_ext, na.rm = TRUE),
                     max(out_dat$spat_ext, na.rm = TRUE),
                     length = 200)
    }
    
    out_pred <- out_dat %>%
      dplyr::group_by(spat_ext) %>%
      dplyr::summarise(samp_eff = median(samp_eff, na.rm = TRUE), .groups = "drop") %>%
      dplyr::right_join(data.frame(spat_ext = spatvec), by = "spat_ext") %>%
      tidyr::fill(samp_eff, .direction = "updown")
  }
  
  pred <- predict(scam1, newdata = out_pred, se = TRUE, type = "response")
  
  out_pred <- out_pred %>%
    dplyr::mutate(S = pred$fit,
                  S_low = pred$fit - 2 * pred$se.fit,
                  S_high = pred$fit + 2 * pred$se.fit)
  
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
  data <- sSBR_object$data
  smooth <- sSBR_object$smooth
  method <- attr(sSBR_object, "method")

  x_dat <- if ("spat_ext" %in% names(data)) "spat_ext" else "effort"
  x_sm <- if ("spat_ext" %in% names(smooth)) "spat_ext" else "effort"

  xlab <- switch(
    method,
    euclidean = "Euclidean distance",
    toroidal = "Toroidal distance",
    area = "Cumulative convex hull area",
    "Spatial extent"
  )

  graphics::plot(
    NA,
    xlim = range(data[[x_dat]], na.rm = TRUE),
    ylim = range(data$S, na.rm = TRUE),
    xlab = xlab,
    ylab = "Cumulative species richness",
    las = 1,
    ...
  )

  if (all_lines) {
    a <- 1 / sqrt(nrow(data) / 300)
    alpha <- pmax(0.05, pmin(0.75, a))

    for (i in unique(data$id)) {
      tmp <- data[data$id == i, , drop = FALSE]
      graphics::lines(
        tmp[[x_dat]],
        tmp$S,
        col = colorspace::adjust_transparency(col, alpha = alpha)
      )
    }
  }

  graphics::polygon(
    c(smooth[[x_sm]], rev(smooth[[x_sm]])),
    c(smooth$S_low, rev(smooth$S_high)),
    col = colorspace::adjust_transparency(col, alpha = 0.5),
    border = NA
  )

  graphics::lines(
    smooth[[x_sm]],
    smooth$S,
    col = colorspace::darken(col, amount = 0.4),
    lwd = 3
  )

  invisible(sSBR_object)
}

