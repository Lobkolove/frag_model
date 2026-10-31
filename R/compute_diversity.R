compute_diversity <- function(
  data,
  steps = c("post_fragmentation", "final"),
  frag_levels = c("Low" = 0.2, "Medium" = 0.5, "High" = 0.8),
  metadata_cols = NULL,
  species_col = "species_id",
  abundance_col = "n",
  sample_col = "sample_id",
  step_col = "step_label",
  sim_col = "sim_id",
  fragmentation_col = "fragmentation",
  step_labels = c(
    "Post-fragmentation",
    "End of simulation"
  )
) {

  # Input validation
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }

  if (length(step_labels) != length(steps)) {
    stop("`step_labels` must have the same length as `steps`.")
  }

  if (is.null(metadata_cols)) {
    metadata_cols <- character(0)
  }

  if (!is.character(metadata_cols)) {
    stop("`metadata_cols` must be NULL or a character vector of column names.")
  }

  required_cols <- c(
    species_col,
    abundance_col,
    sample_col,
    step_col,
    sim_col,
    fragmentation_col,
    metadata_cols
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

  if (nrow(data) == 0) {
    stop("No rows remain after removing missing species identifiers.")
  }

  # Dataset-level metadata must have one unique combination
  metadata <- data |>
    dplyr::select(dplyr::all_of(metadata_cols)) |>
    dplyr::distinct()

  if (nrow(metadata) > 1) {
    stop(
      "The supplied metadata columns do not have one unique combination ",
      "within this dataset. Multiple combinations were found for: ",
      paste(metadata_cols, collapse = ", ")
    )
  }

  # Handle the no-metadata case
  if (length(metadata_cols) == 0) {
    metadata <- tibble::tibble(.metadata_dummy = 1L)
  }

  # Values assumed to be constant within a dataset
  sim_id <- data[[sim_col]][1]
  fragmentation <- data[[fragmentation_col]][1]

  purrr::map2_dfr(
    steps,
    step_labels,
    \(current_step, current_label) {

      subset <- data |>
        dplyr::filter(.data[[step_col]] == current_step)

      if (nrow(subset) == 0) {

        results <- tibble::tibble(
          scale = factor(
            c("sample", "landscape"),
            levels = c("sample", "landscape")
          ),
          richness = NA_real_,
          hill_shannon = NA_real_,
          hill_simpson = NA_real_,
          evenness = NA_real_
        )

      } else {

        # Site-by-species abundance matrix
        spec_table <- subset |>
          dplyr::select(
            dplyr::all_of(
              c(sample_col, species_col, abundance_col)
            )
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

        # Alpha richness and Hill numbers
        sample_richness <- rowSums(spec_table > 0)
        sample_shannon <- as.numeric(vegan::renyi(spec_table, scales = 1, hill = TRUE))
        sample_simpson <- as.numeric(vegan::renyi(spec_table, scales = 2, hill = TRUE))

        alpha <- tibble::tibble(
          scale = factor(
            "sample",
            levels = c("sample", "landscape")
          ),
          richness = mean(sample_richness),
          hill_shannon = mean(sample_shannon),
          hill_simpson = mean(sample_simpson),
          evenness = mean(sample_shannon / sample_richness, na.rm = TRUE)
        )

        # Gamma richness and Hill numbers
        pooled <- colSums(spec_table)
        gamma <- tibble::tibble(
          scale = factor(
            "landscape",
            levels = c("sample", "landscape")
          ),
          richness = sum(pooled > 0),
          hill_shannon = as.numeric(vegan::renyi(pooled, scales = 1, hill = TRUE)),
          hill_simpson = as.numeric(vegan::renyi(pooled, scales = 2, hill = TRUE)), 
          evenness = hill_shannon / richness
        )

        results <- dplyr::bind_rows(alpha, gamma)
      }

      # Add simulation metadata and reshape indices
      long_results <- results |>
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
          ),
          scale = factor(
            scale,
            levels = c("sample", "landscape"),
            labels = c("Local scale", "Landscape scale")
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
          cols = c(
            richness,
            hill_shannon,
            hill_simpson,
            evenness
          ),
          names_to = "index",
          values_to = "value"
        )

      # Attach one-row dataset metadata to every result row
      long_results <- tidyr::crossing(
        long_results,
        metadata
      )

      long_results
    }
  ) |>
    dplyr::select(
      -dplyr::any_of(".metadata_dummy"),
      dplyr::all_of(metadata_cols),
      dplyr::everything()
    )
}
