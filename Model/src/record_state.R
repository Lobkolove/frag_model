record_step <- function(
  grid,
  agents,
  agents_grid,
  step
) {
  structure(
    list(
      grid = grid,
      agents = agents,
      agents_grid = agents_grid,
      step = step
    ),
    class = c("model_state", "core_state", "list")
  )
}


# Add metadata and species abundance table to an existing core_state.
record_state <- function(core_state,
                         sim_id,
                         master_seed,
                         grid_size,
                         ac_amount,
                         fragmentation = NA_real_,
                         habitat = NA_real_,
                         step_label = NULL) {
  
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
      sim_id = sim_id,
      master_seed = master_seed,
      step = core_state$step,
      step_label = step_label,
      ac_amount = ac_amount,
      fragmentation = fragmentation,
      habitat = habitat,
      grid_size = grid_size,
      grid = core_state$grid,
      agents = core_state$agents,
      agents_grid = core_state$agents_grid,
      ss_abund = ss_abund
    ),
    class = c("model_state", "full_state", "list")
  )
}
