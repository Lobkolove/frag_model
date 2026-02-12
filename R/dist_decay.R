#' Distance-decay of similarity
#'
#' Estimate pairwise similarity of communities as a function of spatial distance.
#'
#' @param model_sample A data frame of model sample output with columns:
#'                     - species abundances (prefixed with "sp")  
#'                     - coordinates (x_loc, y_loc)
#' @param binary Logical; if TRUE, abundance data are converted to presence/absence before computing similarity.
#' @param method Character; (dis)similarity index to use (see \code{\link[vegan]{vegdist}}).
#'
#' @return A list with components:
#'   - \code{data}: raw pairwise distances and similarities
#'   - \code{smooth}: GAM-smoothed similarity with confidence intervals
#' @export
#' 
#' @noRd
dist_decay <- function(model_sample, 
                       binary = FALSE, 
                       dist_type = c("euclidean", "toroidal"),
                       method = "bray") { 
  
  dist_type <- match.arg(dist_type)

  # Drop empty sites
  model_sample <- model_sample %>%
    dplyr::filter(rowSums(across(starts_with("sp"))) > 0)
  
  # Extract community matrix and drop empty sites
  comm <- model_sample %>%
    dplyr::select(starts_with("sp"))
  
  # Generate matrix with pairwise spatial distance between sites 
  coords <- model_sample %>% 
    dplyr::select(x_loc, y_loc) 
  
  if (dist_type == "euclidean") {
    d <- stats::dist(coords) # spatial Euclidean distances
  } else if (dist_type == "toroidal") {
    d <- toroidal_dist(coords, model_sample$grid_size[1]) # toroidal distances 
  }
  
  similarity <- 1 - vegan::vegdist(comm, 
                                   method = method,
                                   binary = binary) 
  similarity[!is.finite(similarity)] <- NA
  
  out_dat <- data.frame(distance = as.numeric(d),
                        similarity = as.numeric(similarity))
  
  # order by increasing distance
  out_dat <- out_dat[order(out_dat$distance), ]
  
  # Toroidal cutoff: distances beyond half the grid size are not meaningful
  if (dist_type == "euclidean") {
    d_cut <- model_sample$grid_size[1] / 2
    out_dat <- out_dat %>%
      dplyr::filter(distance <= d_cut)
  }
  
  out_pred <- data.frame(distance   = seq(min(out_dat$distance),
                                          max(out_dat$distance),
                                          length = 200),
                         similarity = NA)
  
  # Fit model - GAM with monotonously increasing constraint
  gam1 <- mgcv::gam(similarity ~ s(distance), data = out_dat)
  
  # Predictions - SCAM
  pred <- stats::predict(gam1, out_pred, se = T)
  
  out_pred$similarity <- pred$fit
  out_pred$simi_low   <- pred$fit - 2*pred$se.fit
  out_pred$simi_high  <- pred$fit + 2*pred$se.fit
  
  out <- list(data   = out_dat,    # distance, similarity
              smooth = out_pred)    # distance, similarity, CI
  
  class(out) <- "dist_decay"
  attr(out, "dissimilarity_method") <- method 
  attr(out, "distance_type") <- dist_type
  return(out)
  
}

#' @export
plot.dist_decay <- function(dd_object,
                            col = "violetred4",
                            ...) {
  
  dat <- dd_object$data
  sm <- dd_object$smooth
  
  method <- attr(dd_object, "dissimilarity_method")
  dist_type <- attr(dd_object, "distance_type")
  
  # Set transparency level based on number of observations
  a <- 1 / sqrt(nrow(dat) / 300)
  alpha <- pmax(0.05, pmin(1, a))
  
  # Scatterplot
  graphics::plot(dat$distance,
                 dat$similarity,
                 col = colorspace::adjust_transparency(col, alpha = alpha),
                 pch = 16,
                 cex = .75,
                 xlab = paste0(dist_type, " distance"),
                 ylab = paste0("Similarity (1 - ", method, " dissimilarity)"))
  
  
  # Confidence ribbon
  graphics::polygon(
    c(sm$distance, rev(sm$distance)),
    c(sm$simi_low, rev(sm$simi_high)),
    col = colorspace::adjust_transparency(col, alpha = 0.25),
    border = NA
  )
  
  # Prediction line
  graphics::lines(sm$distance,
                  sm$similarity,
                  col = colorspace::darken(col, amount = 0.4),
                  lwd = 3)
  
  invisible(dd_object)
}

