################ Simulation initialization Function ############################

# The function combines 3 functions ("generate_grid", "generate_agents", and 
# "distribute_agent") together to simplify the simulation initialization process.
# The function then returns the simulation grid, the updated agents list and 
# the agents distribution grid in the form of a list

initialize <- function(
  grid_size,
  n_species,
  n_pop,
  ac_amount, 
  niche_breadth,
  master_seed = NULL,
  seed_landscape = NULL,
  seed_distribution = NULL,
  random_distribution = FALSE,
  track_ad = FALSE
) {
  
  grid <- fbm_fft(gr_size = grid_size,
                  ac_amount = ac_amount,
                  seed = seed_landscape) 
  
  agents_init <- generate_agents(n_species = n_species,
                                 pop = n_pop, 
                                 grid = grid,
                                 seed = seed_distribution) 
  
  agents_grid_and_list <- distribute_agent(gr_size = grid_size,
                                           agents_list = agents_init, 
                                           space = grid, 
                                           nb = niche_breadth,
                                           random_distribution = random_distribution) 
  
  agents_grid <- if (track_ad) agents_grid_and_list$agents_grid else NULL
  agents <- agents_grid_and_list$agents_list 
  
  state_out <- record_step(
    grid = grid,
    agents = agents,
    agents_grid = agents_grid,
    step = 0L
  )
  return(state_out)
}