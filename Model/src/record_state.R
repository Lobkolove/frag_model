record_step <- function(
  grid,
  agents,
  step
) {
  structure(
    list(
      grid = grid,
      agents = agents,
      step = step
    ),
    class = c("model_state", "core_state", "list")
  )
}


# Add metadata and species abundance table to an existing core_state.
record_state <- function(core_state,
                         meta = list(),
                         step_label) {
  
  ss_abund <- if (nrow(core_state$agents) > 0) {
    core_state$agents |>
      dplyr::count(x_loc, y_loc, species_id, name = "n")
  } else {
    tibble::tibble(
      x_loc = integer(),
      y_loc = integer(),
      species_id = integer(),
      n = integer()
    )
  }

  structure(
    list(
      meta = meta,
      step = core_state$step,
      step_label = step_label,
      grid = core_state$grid,
      agents = core_state$agents,
      ss_abund = ss_abund
    ),
    class = c("model_state", "full_state", "list")
  )
}
