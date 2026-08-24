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
dist_decay <- function(
  data,
  input_format = c("long", "wide"),
  dist_type = c("euclidean", "toroidal"),
  distvec = NULL,
  cutoff = NULL,
  method = "bray",
  binary = FALSE,
  CI = FALSE
) { 
  
  # Input validation
  input_format <- match.arg(input_format)
  dist_type <- match.arg(dist_type)

  # If CIs are requested, warn that they are likely inappropriate due to dependence between points
  if (CI) {
    warning("Confidence intervals are likely inappropriate due to dependence between points.")
  }

  if (method == "bray" && binary) {
    message("Note that with binary = TRUE, the Bray-Curtis dissimilarity is equivalent to the Sørensen index.")
  } else if (method == "jaccard" && !binary) {
    message("Note that with binary = FALSE, the Jaccard dissimilarity is equivalent to the Ruzicka index.")
  }

  # Access grid size
  grid_size <- unique(data$grid_size)

  # If input_format is "long", convert to wide format
  if (input_format == "long") {
    data <- data %>%
      dplyr::filter(!is.na(species_id)) %>%
      tidyr::pivot_wider(
        id_cols = c("sim_id", "step_label", "fragmentation", "x_loc", "y_loc"),
        names_from = "species_id",
        values_from = "n",
        values_fill = 0,
        names_prefix = "sp_"
      )
    } else {
      # Drop empty sites
      data <- data %>%
        dplyr::filter(rowSums(across(starts_with("sp"))) > 0)
      
      # If sp_NA column exists, drop it
      if ("sp_NA" %in% colnames(data)) {
        data <- data %>%
          dplyr::select(-sp_NA)
      }
  }

  # Extract community matrix and drop empty sites
  comm <- data %>%
    dplyr::select(starts_with("sp"))

  # Generate matrix with pairwise spatial distance between sites 
  coords <- data %>% 
    dplyr::select(x_loc, y_loc) 
  
  if (dist_type == "euclidean") {
    d <- stats::dist(coords) # spatial Euclidean distances
  } else if (dist_type == "toroidal") {
    d <- toroidal_dist(coords, grid_size) # toroidal distances 
  }
  
  # Compute pairwise similarity (1 - dissimilarity) between sites
  similarity <- 1 - vegan::vegdist(comm, method = method, binary = binary) 
  similarity[!is.finite(similarity)] <- NA
  
  # Create data frame with distances and similarities, ordered by increasing distance
  out_dat <- data.frame(distance = as.numeric(d),
                        similarity = as.numeric(similarity))
  out_dat <- out_dat[order(out_dat$distance), ]
  
  # Toroidal cutoff: distances beyond half the grid size are not meaningful
  if (is.null(cutoff)) {
    if (!is.null(distvec)) cutoff <- max(distvec) else cutoff <- grid_size / 2
  }
  if (cutoff > grid_size / 2) {
    warning("Cutoff distance exceeds half the grid size, which might lead to unreliable results.")
  }
  out_dat <- out_dat %>%
    dplyr::filter(distance <= cutoff)

  # If distvec is not provided, create a default sequence of distances for predictions
  if (is.null(distvec)) {
    distvec <- seq(
      from = min(out_dat$distance),
      to = max(out_dat$distance),
      length = 200
    )
  }
  
  # Create a data frame for predictions with the specified distance vector
  out_pred <- tibble(distance = distvec)
  
  # Fit model (GAM) to the data
  gam1 <- mgcv::gam(similarity ~ s(distance), data = out_dat)
  
  # Predictions
  pred <- stats::predict(gam1, out_pred, se = TRUE)
  
  # Add predictions (and confidence intervals) to prediction tibble
  out_pred$similarity <- pred$fit
  if (CI) {
    out_pred$simi_low   <- pred$fit - 2*pred$se.fit
    out_pred$simi_high  <- pred$fit + 2*pred$se.fit
    # These confidence intervals ignore the dependence between the points,
    # so they are likely inappropiate.
  } 
  
  # Return a list with raw data and smoothed predictions
  out <- list(data   = out_dat,    # distance, similarity
              smooth = out_pred)    # distance, similarity, CI
  
  class(out) <- "dist_decay"
  attr(out, "dissimilarity_method") <- method 
  attr(out, "distance_type") <- dist_type
  return(out)
  
}

grouped_ddecay <- function(
  model_sample,
  input_format = c("long", "wide"),
  group_cols = c("fragmentation", "step_label"),
  dist_type = c("euclidean", "toroidal"),
  distvec = NULL,
  method = "bray",
  binary = FALSE,
  on_error = c("stop", "warn")
) {

  dist_type <- match.arg(dist_type)
  on_error <- match.arg(on_error)

  if (is.null(distvec)) {
    message("Note that for grouped distance decay, providing a custom `distvec` is recommended to ensure consistent distance values across groups.")
  }

  if (!is.data.frame(model_sample)) {
    stop("`model_sample` must be a data frame.")
  }

  if (!is.character(group_cols) || length(group_cols) == 0) {
    stop("`group_cols` must contain at least one column name.\nIf you want to compute distance decay for the entire dataset, use `dist_decay()` instead.")
  }

  missing_cols <- setdiff(group_cols, names(model_sample))

  if (length(missing_cols) > 0) {
    stop(
      "The following grouping columns are missing from the input data: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  model_sample |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(
      \(data_group, group_keys) {
        result <- tryCatch(
          {
            decay <- dist_decay(
              data = data_group,
              input_format = input_format,
              dist_type = dist_type,
              distvec = distvec,
              method = method,
              binary = binary
            )

            predictions <- decay$smooth |>
              dplyr::select(distance, similarity)

            predictions
          },
          error = function(e) {
            if (on_error == "stop") {
              stop(e)
            }

            warning(
              "Distance decay failed for group: ",
              paste(
                paste0(names(group_keys), "=", unlist(group_keys)),
                collapse = ", "
              ),
              "\nReason: ",
              conditionMessage(e)
            )

            tibble::tibble(
              distance = distvec,
              similarity = NA_real_
            )
          }
        )

        result
      },
      .keep = TRUE
    ) |>
    dplyr::ungroup()

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