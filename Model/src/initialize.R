################ Simulation initialization Function ############################

# The function combines 3 functions ("generate_grid", "generate_agents", and 
# "distribute_agent") together to simplify the simulation initialization process.
# The function then returns the simulation grid, the updated agents list and 
# the agents distribution grid in the form of a list

initialize <- function(frag, hab, ac, nb) {
  
  grid <- generate_grid(gr_size = mod_par$grid_size, 
                        ac_amount = ac, 
                        frag_amount = frag,
                        hab_amount =  hab) 
  
  agents_init <- generate_agents(n_species = mod_par$n_species,
                                 pop = mod_par$n_pop, grid) 
  
  agents_grid_and_list <- distribute_agent(gr_size = mod_par$grid_size,
                                           agents_list = agents_init, 
                                           space = grid, 
                                           nb = nb) 
  
  agents_grid <- agents_grid_and_list$agents_raster
  agents <- agents_grid_and_list$agents_list 
  
  
  return_list <- list(grid = grid, agents = agents, agents_grid = agents_grid)
  return(return_list)
}