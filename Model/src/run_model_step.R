# Core model logic for a single timestep
run_model_step <- function(model_state,
                           var_par,
                           switch) {
  
  grid <- model_state$grid
  agents <- model_state$agents
  agents_grid <- model_state$agents_grid

  # Birth
  step1 <- birth(
    agents = agents,
    agents_grid = agents_grid,
    grid = grid,
    nb = var_par$nb,
    disp = var_par$disp,
    d_dis = var_par$disp_dist
  )
  agents <- step1$agents
  agents_grid <- step1$agents_grid
  
  # Death
  step2 <- death(
    agents = agents,
    agents_grid = agents_grid,
    grid = grid,
    edge_fac = var_par$edge
  )
  agents <- step2$agents
  agents_grid <- step2$agents_grid
  
  # Immigration
  if (switch$immigration == 1) {
    step3 <- immigration(
      agents = agents, 
      agents_grid = agents_grid, 
      grid = grid)
    agents <- step3$agents
    agents_grid <- step3$agents_grid
  }
  
  step_out <- record_step(grid = grid,
                          agents = agents,
                          agents_grid = agents_grid,
                          step = model_state$step + 1L)
  
  if ("core_state" %in% class(model_state)) {
    return(step_out)
  } else if ("full_state" %in% class(model_state)) {
    full_state <- record_state(
      core_state = step_out,
      sim_id = model_state$sim_id,
      master_seed = model_state$master_seed,
      grid_size = model_state$grid_size,
      ac_amount = model_state$ac_amount,
      fragmentation = model_state$fragmentation,
      habitat = model_state$habitat,
      step_label = NULL
    )
    return(full_state)
  }
}  