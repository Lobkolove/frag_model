################## Dynamic model function reduced to core ######################

# Reduced version of GeDo_run.R which only runs the core model logic.
# It returns complete system snapshot at selected timesteps, in form of a list 
# which includes the following:
# landscape grid (grid), agents list (agents) and a data frame with aggregated
# species abundances for each site (ss_abund).

clean_run <- function(mod_par,
                      var_par,
                      switch,
                      sim_id,
                      record_steps = c("start", 
                                       "pre_fragmentation", 
                                       "post_fragmentation", 
                                       "final"),
                      seed = NULL) {
  
  force(sim_id)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Initialization ----------------------------------------------------------
  
  species_sequence <- 1:mod_par$n_species
  steps_1 <- mod_par$steps_pre_frag
  steps_2 <- mod_par$steps_post_frag
  
  
  # Initialize model with 100% habitat
  model_start <- initialize(
    ac   = var_par$ac,
    nb   = var_par$nb
  )
  
  # extract simulation space, agents grid and agents list
  grid         <- model_start$grid
  agents_grid  <- model_start$agents_grid
  agents       <- model_start$agents
  
  # Storage for recorded states
  state_list <- list()
  
  # Helper to store a full system snapshot
  record_state <- function(step_label, step_number, fragmented = TRUE) {
    
    # Compute abundances per site per species (makes later sampling much faster)
    ss_abund <- agents %>%
      dplyr::count(x_loc, y_loc, species_id, name = "n")
    
    list(
      sim_id        = sim_id,
      master_seed   = master_seed,
      step          = step_number,
      step_label    = step_label,
      fragmentation = ifelse(fragmented, var_par$frag, NA),
      ac_amount     = var_par$ac,
      habitat       = ifelse(fragmented, var_par$hab, NA),
      grid_size     = mod_par$grid_size,
      grid          = grid,
      agents        = agents,
      agents_grid   = agents_grid,
      ss_abund      = ss_abund
    )
  }
  
  # Record before the first step
  if ("start" %in% record_steps) {
    state_list[["start"]] <- record_state(
      step_label  = "start",
      step_number = 0,
      fragmented = FALSE
    )
  }

  # Pre-fragmentation -------------------------------------------------------
  
  for (i in seq_len(steps_1)) {
    
    start.time <- Sys.time()
    
    if (nrow(agents) > 0) {
      
      # Run a full model step and update agents list and grid
      step_out <- run_model_step(
        grid         = grid,
        agents       = agents,
        agents_grid  = agents_grid,
        var_par      = var_par,
        switch       = switch
      )
      agents <- step_out$agents
      agents_grid <- step_out$agents_grid
      
    } else {
      message("ALL DEAD before fragmentation")
      break
    }
    
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time, 2)
    if (switch$print_agents == 1) {
      cat("step", i, "took", time.taken, "with", nrow(agents), "agents\n")
    }
  }
  
  # Record right before fragmentation
  if ("pre_fragmentation" %in% record_steps) {
    state_list[["pre_fragmentation"]] <- record_state(
      step_label  = "pre_fragmentation",
      step_number = steps_1,
      fragmented = FALSE
    )
  }

  # Fragmentation event -----------------------------------------------------

  current_state <- record_state(
    step_label  = "pre_fragmentation",
    step_number = steps_1
  )

  frag_out <- fragment(full_state = current_state,
                       habitat = var_par$hab,
                       fragmentation = var_par$frag)
  
  grid        <- frag_out$grid
  agents      <- frag_out$agents
  agents_grid <- frag_out$agents_grid
  
  if (switch$print_agents == 1) {
    print("FRAGMENTATION")
  }
  
  # Record immediately after fragmentation
  if ("post_fragmentation" %in% record_steps) {
    state_list[["post_fragmentation"]] <- frag_out
  }
  

  # Post-fragmentation ------------------------------------------------------

  for (j in seq_len(steps_2)) {
    
    if (nrow(agents) > 0) {
      
      # Run a full model step and update agents list and grid
      step_out <- run_model_step(
        agents       = agents,
        agents_grid  = agents_grid,
        grid         = grid,
        var_par      = var_par,
        switch       = switch
      )
      agents <- step_out$agents
      agents_grid <- step_out$agents_grid
      
    } else {
      message("ALL DEAD after fragmentation")
      break
    }
  }
  
  # Record final state
  if ("final" %in% record_steps) {
    state_list[["final"]] <- record_state(
      step_label  = "final",
      step_number = steps_1 + steps_2
    )
  }
  
  return(state_list)
}
