spat_acc <- function(
  sampled_data,
  method = c("area", "euclidean", "toroidal"),
  n_focal = NULL,
  compute_richness = FALSE
) {

  method <- match.arg(method)

  # numeric convex hull area via shoelace formula (no sf)
  hull_area_numeric <- function(xy) {
    if (nrow(xy) < 3) {
      return(0)
    }
    h <- grDevices::chull(xy[, 1], xy[, 2])
    h <- c(h, h[1])
    x <- xy[h, 1]
    y <- xy[h, 2]
    abs(sum(x[-1] * y[-length(y)] - x[-length(x)] * y[-1])) / 2
  }

  n <- nrow(sampled_data)
  coords <- sampled_data |>
    dplyr::select(x_loc, y_loc)

  if (compute_richness) {
    sp_cols <- names(sampled_data)[startsWith(names(sampled_data), "sp")]
    if (length(sp_cols) == 0) {
      stop("No species columns found. Expected columns prefixed with 'sp'.")
    }

    pa_table <- sampled_data |>
      dplyr::select(dplyr::all_of(sp_cols)) |>
      dplyr::select(dplyr::where(~ sum(.x, na.rm = TRUE) > 0)) |>
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.integer(.x > 0)))
  }

  pair_dist <- switch(
    method,
    euclidean = as.matrix(stats::dist(coords)),
    toroidal = toroidal_dist(coords, sampled_data$grid_size[1]),
    area = as.matrix(stats::dist(coords))
  )

  focal_ids <- seq_len(n)
  if (!is.null(n_focal) && is.finite(n_focal) && n_focal < n) {
    focal_ids <- sample.int(n, n_focal)
  }

  out_list <- vector("list", length(focal_ids))

  for (k in seq_along(focal_ids)) {
    i <- focal_ids[k]

    dist_to_site <- pair_dist[i, ]

    # random tie-break among equal distances
    new_order <- sample.int(n)
    new_order <- new_order[order(dist_to_site[new_order])]
    new_order <- c(i, new_order[new_order != i])

    if (method %in% c("euclidean", "toroidal")) {
      spat_ext <- dist_to_site[new_order]
    } else {
      coords_ordered <- as.matrix(coords[new_order, , drop = FALSE])

      spat_ext <- numeric(n)
      spat_ext[1] <- 0

      # as in mobr
      spat_ext[2] <- dist_to_site[new_order][2]

      for (j in 3:n) {
        spat_ext[j] <- hull_area_numeric(coords_ordered[1:j, , drop = FALSE])
      }
    }

    if (compute_richness) {
      comm_ordered <- pa_table[new_order, , drop = FALSE]
      comm_bool <- (comm_ordered == 0) * 1
      rich <- apply(comm_bool, 2, cumprod)
      S <- rowSums(1 - rich)

      out_list[[k]] <- tibble::tibble(
        id = i,
        spat_ext = spat_ext,
        samp_eff = seq_len(n),
        S = S
      )
    } else {
      out_list[[k]] <- tibble::tibble(
        sim_id = sampled_data$sim_id[1],
        spat_ext = spat_ext,
        samp_eff = seq_len(n)
      )
    }
  }

  dplyr::bind_rows(out_list)
}

ref_curve <- function(
  data_list,
  quantile = 0.5,
  spatpar = c("area", "euclidean", "toroidal"),
  spatvec = NULL,
  n_focal = NULL
) {

  spatpar <- match.arg(spatpar)

  # Process each tibble separately (safe memory usage)
  out_dats <- map(
    data_list,
    ~ spat_acc(.x, method = spatpar, n_focal = n_focal, compute_richness = FALSE)
  )

  # Pool results (tiny memory footprint)
  pooled_out <- bind_rows(out_dats)

  # Build single shared reference curve
  eff_ref <- pooled_out %>%
    group_by(spat_ext) %>%
    summarise(
      samp_eff = quantile(samp_eff, probs = quantile, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(samp_eff))

  # Interpolate to spatvec if provided
  if (!is.null(spatvec)) {
    eff_ref <- eff_ref %>%
      right_join(tibble(spat_ext = spatvec), by = "spat_ext") %>%
      tidyr::fill(samp_eff, .direction = "updown")
  }

  # Metadata
  attr(eff_ref, "quantile") <- quantile
  class(eff_ref) <- c("effort_reference", class(eff_ref))

  return(eff_ref)
}