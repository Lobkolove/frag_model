#' Distance-decay of similarity
#'
#' Estimate pairwise similarity of communities as a function of spatial distance.
#'
#' @param model_sample A data frame of model sample output with columns:
#'                     - species abundances (prefixed with "sp")  
#'                     - coordinates (loc_x, loc_y)
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
                       method = "bray") { 
  
  # Drop empty sites
  model_sample <- model_sample %>%
    filter(rowSums(across(starts_with("sp"))) > 0)
  
  # Extract community matrix and drop empty sites
  comm <- model_sample %>%
    select(starts_with("sp"))
  
  # Generate matrix with pairwise spatial distance between sites 
  coords <- model_sample %>% 
    select(loc_x, loc_y) 
  
  d <- stats::dist(coords) # spatial Euclidean distance
  
  similarity <- 1 - vegan::vegdist(comm, 
                                   method = method,
                                   binary = binary) 
  similarity[!is.finite(similarity)] <- NA
  
  out_dat <- data.frame(distance = as.numeric(d),
                        similarity = as.numeric(similarity))
  
  # order by increasing distance
  out_dat <- out_dat[order(out_dat$distance), ]
  
  out_pred <- data.frame(distance   = seq(min(out_dat$distance),
                                          max(out_dat$distance),
                                          length = 200),
                         similarity = NA)
  
  # Fit model - GAM with monotonously increasing constraint
  gam1 <- mgcv::gam(similarity ~ s(distance), data = out_dat)
  
  # Predictions - SCAM
  pred <- stats::predict(gam1, out_pred, se = T)
  
  out_pred$similarity <- pred$fit
  out_pred$simi_low <- pred$fit - 2*pred$se.fit
  out_pred$simi_high <- pred$fit + 2*pred$se.fit
  
  out <- list(data   = out_dat,    # distance, similarity
              smooth = out_pred)    # distance, similarity, CI
  
  class(out) <- "dist_decay"
  return(out)
  
}

#' @export
plot.dist_decay <- function(dd_object,
                            col = "red",
                            ...) {
  
  dat <- dd_object$data
  sm <- dd_object$smooth
  
  # Scatterplot
  graphics::plot(dat$distance,
                 dat$similarity,
                 col = adjustcolor(col, alpha.f = .5),
                 pch = 16,
                 cex = .75,
                 xlab = "Spatial distance",
                 ylab = "Similarity",
                 main = "Distance decay")
  
  
  # Confidence ribbon
  graphics::polygon(
    c(sm$distance, rev(sm$distance)),
    c(sm$simi_low, rev(sm$simi_high)),
    col = adjustcolor(col, alpha.f = .25),
    border = NA
  )
  
  # Prediction line
  graphics::lines(sm$distance,
                  sm$similarity,
                  col = col,
                  lwd = 3)
  
  invisible(dd_object)
}

