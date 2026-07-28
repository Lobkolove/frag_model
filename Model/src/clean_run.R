################## Dynamic model function reduced to core ######################

# Reduced version of GeDo_run.R which only runs the core model logic.
# It returns complete system snapshot at selected timesteps, in form of a list
# which includes the following:
# landscape grid (grid), agents list (agents) and a data frame with aggregated
# species abundances for each site (ss_abund).

clean_run <- function(
  mod_par,
  var_par,
  species_par,
  switch,
  sim_id,
  record_steps = c("start", "pre_fragmentation", "post_fragmentation", "final"),
  seed = NULL
) {

  force(sim_id)

  # If seed is provided explicitly, use it. Otherwise, check if a master seed was sourced.
  # If neither is provided, generate a random master seed and set it.
  if (!is.null(seed)) {
    set.seed(seed)
    master_seed <- seed
  } else if (exists("master_seed", envir = .GlobalEnv)) {
    set.seed(master_seed)
  } else {
    master_seed <- round(runif(1, 0, 1e6))
    set.seed(master_seed)
  }
  
  # Check if species specific parameters are enabled
  species_specific <- isFALSE(switch$species_specific_par == 0)
  if (species_specific) {
    message("Species-specific parameters are enabled. `var_par` values will be ignored for species-specific parameters.")
  } else {
    if (switch$dispersal_type == 0) {
      message("Random dispersal is enabled. `var_par$disp` and `var_par$disp_dist` will be ignored.")
    }
  }

  # Check which dispersal type is selected
  dispersal_type <- ifelse(switch$dispersal_type == 0, "random", "short_long")
  dispersal_kernel <- ifelse(switch$kernel_type == 0, "log-normal", "exponential")

  # Initialize metadata for recording
  # (habitat and fragmentation will be updated after fragmentation event)
  meta <- build_meta(
    sim_id = sim_id,
    run_date = Sys.Date(),
    master_seed = master_seed,
    grid_size = mod_par$grid_size,
    ac_amount = var_par$ac,
    habitat = 1,
    fragmentation = NA_real_,
    niche_breadth = var_par$nb,
    edge_effect = var_par$edge,
    dispersal_type = dispersal_type,
    dispersal_kernel = dispersal_kernel,
    dispersal_ratio = var_par$disp,
    dispersal_dist = var_par$disp_dist
  )


  # Initialization ----------------------------------------------------------

  # Number of steps before and after fragmentation
  steps_1 <- mod_par$steps_pre_frag
  steps_2 <- mod_par$steps_post_frag

  # Set seeds for different processes based on switches
  seed_landscape <- if (switch$random_landscape == 1) master_seed + 1 else NULL
  seed_distribution <- if (switch$random_community == 1) master_seed + 2 else NULL
  seed_fragment <- if (switch$random_post_frag == 1) master_seed + 3 else NULL

  # Initialize model with 100% habitat
  model_start <- initialize(
    grid_size = mod_par$grid_size,
    ac_amount = var_par$ac,
    n_species = mod_par$n_species,
    n_pop = mod_par$n_pop,
    species_par = species_par,
    k_inter = mod_par$k_inter,
    k_intra = mod_par$k_intra,
    niche_breadth = var_par$nb,
    species_specific = species_specific,
    master_seed = master_seed,
    seed_landscape = seed_landscape,
    seed_distribution = seed_distribution
  )

  # Storage for recorded states
  state_list <- list()

  # Record before the first step
  if ("start" %in% record_steps) {
    state_list[["start"]] <- record_state(
      model_state = model_start,
      meta = meta,
      step_label = "start"
    )
  }

  # Pre-fragmentation Loop ---------------------------------------------------

  model_state <- model_start

  for (i in seq_len(steps_1)) {
    start.time <- Sys.time()

    if (nrow(model_state$agents) > 0) {
      # Run a full model step and update state
      step_out <- run_model_step(
        model_state = model_state,
        mod_par = mod_par,
        var_par = var_par,
        species_par = species_par,
        switch = switch
      )
      model_state <- step_out
    } else {
      message("ALL DEAD before fragmentation")
      break
    }

    end.time <- Sys.time()
    time.taken <- round(end.time - start.time, 2)
    if (switch$print_agents == 1) {
      cat("step", model_state$step, "took", time.taken, "with", nrow(model_state$agents), "agents\n")
    }
  }

  full_state <- record_state(
    model_state = model_state,
    meta = meta,
    step_label = "pre_fragmentation"
  )

  # Record right before fragmentation
  if ("pre_fragmentation" %in% record_steps) {
    state_list[["pre_fragmentation"]] <- full_state
  }

  # Fragmentation event -----------------------------------------------------

  frag_out <- fragment(
    full_state = full_state,
    habitat = var_par$hab,
    fragmentation = var_par$frag,
    seed = seed_fragment
  )

  if (switch$print_agents == 1) {
    print("FRAGMENTATION")
  }

  # Update metadata with fragmentation parameters
  meta$fragmentation <- var_par$frag
  meta$habitat <- var_par$hab

  # Record immediately after fragmentation
  if ("post_fragmentation" %in% record_steps) {
    state_list[["post_fragmentation"]] <- frag_out
  }

  # Post-fragmentation ------------------------------------------------------

  model_state <- frag_out

  for (j in seq_len(steps_2)) {
    start.time <- Sys.time()

    if (nrow(model_state$agents) > 0) {
      # Run a full model step and update agents list and grid
      step_out <- run_model_step(
        model_state = model_state,
        mod_par = mod_par,
        var_par = var_par,
        species_par = species_par,
        switch = switch
      )
      model_state <- step_out
    } else {
      message("ALL DEAD after fragmentation")
      break
    }

    end.time <- Sys.time()
    time.taken <- round(end.time - start.time, 2)
    if (switch$print_agents == 1) {
      cat("step", model_state$step, "took", time.taken, "with", nrow(model_state$agents), "agents\n")
    }
  }

  # Record final state
  if ("final" %in% record_steps) {
    state_list[["final"]] <- record_state(
      model_state = model_state,
      meta = meta,
      step_label = "final"
    )
  }

  return(state_list)
}
