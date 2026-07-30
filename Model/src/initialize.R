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
  seed_distribution = NULL,
  state_type = c("core_state", "full_state")
) {

  state_type <- match.arg(state_type)

  if (!is.null(master_seed)) {
    if (is.null(seed_landscape)) {
      seed_landscape <- master_seed + 1
    }
    if (is.null(seed_distribution)) {
      seed_distribution <- master_seed + 2
    }
  }

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

  if (state_type == "core_state") {
    return(state_out)
  } else if (state_type == "full_state") {
    meta <- build_meta(
      sim_id = NA,
      run_date = Sys.Date(),
      master_seed = master_seed,
      grid_size = grid_size,
      ac_amount = ac_amount,
      habitat = 1,
      fragmentation = NA,
      niche_breadth = niche_breadth,
      edge_effect = NA,
      dispersal_type = NA,
      dispersal_kernel = NA,
      dispersal_ratio = NA,
      dispersal_dist = NA
    )
    full_state <- record_state(
      model_state = state_out,
      meta = meta,
      step_label = "start"
    )
  }
}
