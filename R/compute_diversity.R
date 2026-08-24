compute_diversity <- function(
  data,
  steps = c("post_fragmentation", "final"),
  frag_levels = c("low" = 0.2, "medium" = 0.5, "high" = 0.8),
  species_col = "species_id",
  abundance_col = "n",
  sample_col = "sample_id",
  step_col = "step_label",
  sim_col = "sim_id",
  fragmentation_col = "fragmentation",
  step_labels = c("Post-fragmentation", "End of simulation")
) {

  # Check arguments
  stopifnot(is.data.frame(data))

  required_cols <- c(
    species_col,
    abundance_col,
    sample_col,
    step_col,
    sim_col,
    fragmentation_col
  )
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Remove rows without a species identity
  data <- data |>
    dplyr::filter(!is.na(.data[[species_col]]))

  # Metadata shared by all steps, assuming it is constant within a dataset
  sim_id <- data[[sim_col]][1]
  fragmentation <- data[[fragmentation_col]][length(data[[fragmentation_col]])]

  # For each step, compute diversity indices and return a tidy data frame
  purrr::map2_dfr(
    steps,
    step_labels,
    \(current_step, current_label) {
      
      subset <- data |>
        dplyr::filter(.data[[step_col]] == current_step)

      if (nrow(subset) == 0) {
        stop(
          "No rows found for step '",
          current_step,
          "'. Aborting."
        )
      }

      # Landscape-level richness
      gamma <- tibble::tibble(
        scale = factor(
          "landscape",
          levels = c("sample", "landscape"),
          labels = c("Sample scale", "Landscape scale")
        ),
        richness = dplyr::n_distinct(subset[[species_col]])
      )

      # Mean sample-level richness
      alpha <- subset |>
        dplyr::group_by(.data[[sample_col]]) |>
        dplyr::summarise(
          richness = dplyr::n_distinct(.data[[species_col]]),
          .groups = "drop"
        ) |>
        dplyr::summarise(
          richness = mean(richness, na.rm = TRUE)
        ) |>
        dplyr::mutate(
          scale = factor(
            "sample",
            levels = c("sample", "landscape"),
            labels = c("Sample scale", "Landscape scale")
          )
        ) |>
        dplyr::select(scale, richness)

      # Site-by-species abundance matrix
      spec_table <- subset |>
        dplyr::select(
          dplyr::all_of(c(sample_col, species_col, abundance_col))
        ) |>
        tidyr::pivot_wider(
          id_cols = dplyr::all_of(sample_col),
          names_from = dplyr::all_of(species_col),
          values_from = dplyr::all_of(abundance_col),
          values_fill = 0,
          names_prefix = "sp_",
          names_sort = TRUE
        ) |>
        dplyr::select(dplyr::starts_with("sp_"))

      # Calculate Hill numbers
      alpha <- alpha |>
        dplyr::mutate(
          hill_shannon = mean(
            as.numeric(
              vegan::renyi(
                spec_table,
                scales = 2,
                hill = TRUE
              )
            )
          ),
          hill_simpson = mean(
            as.numeric(
              vegan::renyi(
                spec_table,
                scales = 3,
                hill = TRUE
              )
            )
          ),
          evenness = hill_shannon / richness
        )

      gamma <- gamma |>
        dplyr::mutate(
          hill_shannon = as.numeric(
            vegan::renyi(
              colSums(spec_table),
              scales = 2,
              hill = TRUE
            )
          ),
          hill_simpson = as.numeric(
            vegan::renyi(
              colSums(spec_table),
              scales = 3,
              hill = TRUE
            )
          ),
          evenness = hill_shannon / richness
        )
      
      # Combine alpha and gamma results, add metadata, and reshape to tidy format
      dplyr::bind_rows(alpha, gamma) |>
        dplyr::mutate(
          sim_id = sim_id,
          fragmentation = factor(
            fragmentation,
            levels = frag_levels,
            labels = names(frag_levels)
          ),
          step_label = factor(
            current_step,
            levels = steps,
            labels = step_labels
          )
        ) |>
        dplyr::select(
          sim_id,
          fragmentation,
          step_label,
          scale,
          richness,
          hill_shannon,
          hill_simpson,
          evenness
        ) |>
        tidyr::pivot_longer(
          cols = richness:evenness,
          names_to = "index",
          values_to = "value"
        ) |> 
        dplyr::mutate(
          index = factor(
            index,
            levels = c("richness", "hill_shannon", "hill_simpson", "evenness"),
            labels = c(
              "Species richness",
              "Hill Shannon (q = 1)",
              "Hill Simpson (q = 2)",
              "Evenness"
            )
          )
        )
    }
  )
}
