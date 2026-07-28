#' Sample habitat cells from an ABM state
#'
#' Extracts a subset of habitat cells from a model state object and returns
#' per-cell species abundances along with spatial coordinates, patch
#' information, and simulation metadata. Cells can be sampled exhaustively,
#' randomly, or using a deterministic chessboard pattern.
#'
#' @param full_state A named list representing the current ABM state.
#'   At minimum, it must contain:
#'   \itemize{
#'     \item \code{grid}: a \code{RasterLayer} defining habitat (non-\code{NA})
#'       and non-habitat (\code{NA}) cells.
#'     \item \code{ss_abund}: a data table with columns
#'       \code{x_loc}, \code{y_loc}, \code{species_id}, and \code{n}, giving
#'       species abundances per cell.
#'     \item \code{sim_id}, \code{master_seed}, \code{step},
#'       \code{grid_size}, \code{fragmentation}: scalar metadata values.
#'   }
#'
#' @param method Sampling method used to select habitat cells. One of:
#'   \describe{
#'     \item{\code{"all"}}{All habitat cells are returned.}
#'     \item{\code{"random"}}{A random subset of habitat cells is returned
#'       (sampling without replacement).}
#'     \item{\code{"chessboard"}}{Every other habitat cell is selected based on
#'       the parity of row and column indices.}
#'   }
#'
#' @param n_samples Integer number of habitat cells to sample when
#'   \code{method = "random"}. Must be provided for random sampling.
#'
#' @param format Output format. Either \code{"long"} (one row per
#'   cell–species combination) or \code{"wide"} (one row per sampled cell, one
#'   column per species).
#'
#' @param seed Optional integer seed for reproducible random sampling.
#'
#' @details
#' Cell coordinates (\code{x_loc}, \code{y_loc}) refer to row and column indices
#' of the underlying raster grid, as returned by
#' \code{raster::rowColFromCell()}.
#'
#' Patch identifiers and patch sizes are computed using
#' \code{toroidal_clump()}, which assigns patch IDs while accounting for the
#' toroidal nature of the grid.
#'
#' In wide format, missing species abundances are filled with zeros.
#'
#' @return A data frame containing sampled cells and associated data.
#'   The exact structure depends on \code{format}:
#'   \itemize{
#'     \item \code{"long"}: one row per sampled cell and species.
#'     \item \code{"wide"}: one row per sampled cell, with species abundances
#'       spread across columns prefixed with \code{"sp_"}.
#'   }
#'
#' @importFrom rlang arg_match
#' @importFrom raster getValues rowColFromCell freq
#' @importFrom dplyr left_join mutate
#' @importFrom tidyr pivot_wider
#'
#' @export
sample_cells <- function(full_state,
                         method = c("all", "random", "checkerboard"),
                         n_samples = NULL,
                         format = c("wide", "long"),
                         seed = NULL) {
  
  method <- rlang::arg_match(method)
  format <- rlang::arg_match(format)
  
  if (!is.null(seed)) set.seed(seed)
  
  # Extract grid and species abundances per cell from state
  grid            <- full_state$grid
  ss_abund        <- full_state$ss_abund
  
  # Extract IDs and coordinates of habitat cells only
  grid_vals <- raster::getValues(grid)
  habitat_cells <- which(!is.na(grid_vals))
  # Convert cell index to row/col
  coords <- raster::rowColFromCell(grid, habitat_cells)
  
  # Select cells to be sampled
  sampled_cells <- switch(
    method,
    
    all = habitat_cells,
    
    random = {
      if (is.null(n_samples)) {
        stop("n_samples must be provided for random sampling")
      }
      sample(habitat_cells, n_samples)    
    },
    
    checkerboard = {
      # Every other cell based on parity of row + column
      keep <- (coords[, 1] + coords[, 2]) %% 2 == 0
      habitat_cells[keep]
    }
  )
  
  # Assess coordinates of sample cells
  sample_coords <- raster::rowColFromCell(grid, sampled_cells)
  
  samples <- tibble::tibble(
    sample_id = seq_len(length(sampled_cells)),
    cell_id = sampled_cells,
    x_loc = sample_coords[, 1],
    y_loc = sample_coords[, 2]
  )
  
  # Assess patch information
  clumped <- toroidal_clump(grid, directions = 4) 
  patch_freq <- raster::freq(clumped, useNA = "no")
  patch_ids <- patch_freq[, 1]
  patch_sizes <- patch_freq[, 2]
  # Create named vector for patches: name = patch ID, value = patch size
  patches <- setNames(patch_sizes, patch_ids)
  # Add patch info to samples df
  samples$patch_id <- clumped[samples$cell_id]
  samples$patch_size <- unname(patches[as.character(samples$patch_id)])
  
  # Add static metadata columns 
  meta <- full_state$meta
  meta_df <- tibble::tibble(
    sim_id           = meta$sim_id,
    master_seed      = meta$master_seed,
    grid_size        = meta$grid_size,
    ac_amount        = meta$ac_amount,
    habitat          = meta$habitat,
    fragmentation    = meta$fragmentation,
    niche_breadth    = meta$niche_breadth,
    edge_effect      = meta$edge_effect,
    dispersal_type   = meta$dispersal_type,
    dispersal_kernel = meta$dispersal_kernel,
    dispersal_ratio  = meta$dispersal_ratio,
    disp_dist        = meta$dispersal_dist,
    step             = full_state$step,
    step_label       = full_state$step_label,
    samp_method      = method
  )

  # Merge metadata with samples df (metadata first, then sample-specific columns)
  meta_samp <- cbind(meta_df, samples)
  
  
  # Merge samples df and species abundances per cell
  out_long <- dplyr::left_join(meta_samp, ss_abund, by = c("x_loc", "y_loc"))
  
  # Return long format
  if (format == "long") return(out_long)

  # If needed, reformat and return wide format
  out_wide <- out_long %>%
    tidyr::pivot_wider(
      names_from = species_id,
      values_from = n,
      values_fill = 0,
      names_prefix = "sp_",
      names_sort = TRUE
    ) %>%
    dplyr::select(-sp_NA)
  
  return(out_wide)

}
