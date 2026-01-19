# Core model logic for a single timestep
run_model_step <- function(agents,
                           agents_grid,
                           grid,
                           var_par,
                           switch) {
  
  # Birth
  step1 <- birth(
    agents = agents,
    agents_grid = agents_grid,
    grid = grid,
    NB = var_par$nb,
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
    step3 <- immigration(agents, agents_grid, grid)
    agents <- step3$agents
    agents_grid <- step3$agents_grid
  }
  
  list(
    agents = agents,
    agents_grid = agents_grid
  )
}
