# Core model logic for a single timestep
run_model_step <- function(
  model_state,
  var_par,
  switch
) {
  grid <- model_state$grid
  agents <- model_state$agents

  # Birth
  agents <- birth(
    agents = agents,
    grid = grid,
    nb = var_par$nb,
    disp = var_par$disp,
    d_dis = var_par$disp_dist
  )

  # Death
  agents <- death(
    agents = agents,
    grid = grid,
    edge_fac = var_par$edge
  )

  # Immigration
  if (switch$immigration == 1) {
    agents <- immigration(
      agents = agents,
      grid = grid
    )
  }

  step_out <- record_step(
    grid = grid,
    agents = agents,
    step = model_state$step + 1L
  )

  if ("core_state" %in% class(model_state)) {
    return(step_out)
  } else if ("full_state" %in% class(model_state)) {
    full_state <- record_state(
      core_state = step_out,
      meta = model_state$meta,
      step_label = NULL
    )
    return(full_state)
  }
}  