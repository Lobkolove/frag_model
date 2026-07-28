# Core model logic for a single timestep
run_model_step <- function(
  model_state,
  mod_par,
  var_par,
  species_par,
  switch
) {

  grid <- model_state$grid
  agents <- model_state$agents
  cells <- model_state$cells

  d_type <- ifelse(switch$dispersal_type == 1, "short_long", "random")
  kernel <- ifelse(switch$kernel_type == 1, "exponential", "log-normal")
  species_specific <- ifelse(switch$species_specific_par == 0, FALSE, TRUE)

  # Birth
  agents <- birth(
    agents = agents,
    grid = grid,
    k_inter = mod_par$k_inter,
    k_intra = mod_par$k_intra,
    species_par = species_par,
    species_specific = species_specific,
    d_type = d_type,
    kernel = kernel,
    nb = ifelse(species_specific, NULL, var_par$nb),
    disp = ifelse(species_specific, NULL, var_par$disp),
    d_dis = var_par$disp_dist,
    habitat_cells = cells$habitat
  )

  edge_effects <- ifelse(switch$edge_effect == 1, TRUE, FALSE)
  # Death
  agents <- death(
    agents = agents,
    grid = grid,
    species_par = species_par,
    edge_effects = edge_effects,
    edge_cells = cells$edges,
    species_specific = species_specific,
    death_rate = ifelse(species_specific, NULL, mod_par$death_rate),
    edge_factor = ifelse(species_specific, NULL, var_par$edge)
  )

  # Immigration
  if (switch$immigration == 1) {
    agents <- immigration(
      agents = agents,
      grid = grid,
      n_immigrants = mod_par$n_immigrants,
      n_species = mod_par$n_species,
      k_inter = mod_par$k_inter,
      k_intra = mod_par$k_intra,
      species_par = species_par,
      species_specific = species_specific,
      nb = ifelse(species_specific, NULL, var_par$nb)
    )
  }

  step_out <- record_step(
    grid = grid,
    agents = agents,
    step = model_state$step + 1L,
    cells = cells
  )

  if ("core_state" %in% class(model_state)) {
    return(step_out)
  } else if ("full_state" %in% class(model_state)) {
    full_state <- record_state(
      model_state = step_out,
      meta = model_state$meta,
      step_label = NULL
    )
    return(full_state)
  }
}  
