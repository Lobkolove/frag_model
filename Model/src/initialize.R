################ Simulation initialization Function ############################

# The function combines 3 functions ("generate_grid", "generate_agents", and
# "distribute_agent") together to simplify the simulation initialization process.
# The function then returns the simulation grid, the updated agents list and
# the agents distribution grid in the form of a list

initialize <- function(
  grid_size,
  ac_amount,
  n_species,
  n_pop,
  species_par,
  k_inter,
  k_intra,
  niche_breadth = NULL,
  species_specific = FALSE,
  master_seed = NULL,
  seed_landscape = NULL,
  seed_distribution = NULL
) {

  grid <- fbm_fft(
    gr_size = grid_size,
    ac_amount = ac_amount,
    seed = seed_landscape
  )

  cells <- list(
    habitat = 1:grid_size^2,
    edges = NULL
  )

  agents_init <- generate_agents(
    n_species = n_species,
    n_pop = n_pop,
    seed = seed_distribution
  )

  agents <- distribute_agent(
    agents = agents_init,
    grid = grid,
    species_par = species_par,
    k_inter = k_inter,
    k_intra = k_intra,
    nb = niche_breadth,
    species_specific = species_specific,
  )

  state_out <- record_step(
    grid = grid,
    agents = agents,
    step = 0L,
    cells = cells
  )

  return(state_out)
}
