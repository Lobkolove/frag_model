################## Dynamic model function reduced to core ######################

# Reduced version of GeDo_run.R which only runs the model logic and returns a
# complete system state (landscape, agents list and agent grid) at selected
# timesteps.

clean_run <- function(mod_par,
                      var_par,
                      switch,
                      sim_id,
                      record_steps = c("start", "post_fragmentation", "final"),
                      seed = NULL) {
  

  # Initialization ----------------------------------------------------------
  
  species_sequence <- 1:mod_par$n_species
  steps_1 <- mod_par$steps_pre_frag
  steps_2 <- mod_par$steps_post_frag
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Initialize model with 100% habitat
  model_start <- initialize(
    frag = 0,
    hab  = 1,
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
  record_state <- function(step_label, step_number) {
    list(
      sim_id        = sim_id,
      step          = step_number,
      step_label    = step_label,
      fragmentation = var_par$frag,
      habitat       = var_par$hab,
      grid_size     = mod_par$grid_size,
      grid          = grid,
      agents        = agents,
      agents_grid   = agents_grid
    )
  }
  
  # Record immediately after fragmentation
  if ("start" %in% record_steps) {
    state_list[["start"]] <- record_state(
      step_label  = "start",
      step_number = 1
    )
  }

  # Pre-fragmentation -------------------------------------------------------
  
  for (i in 1:steps_1) {
    
    start.time <- Sys.time()
    
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
      message("ALL DEAD before fragmentation")
      break
    }
    
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time, 2)
    if (switch$print_agents == 1) {
      cat("step", i, "took", time.taken, "with", nrow(agents), "agents\n")
    }
  }
  

# Fragmentation event -----------------------------------------------------


  
  frag_out <- cookie_cutting(
    grid        = grid,
    agents      = agents,
    agents_grid = agents_grid,
    hab         = var_par$hab,
    frag        = var_par$frag
  )
  
  grid        <- frag_out$grid
  agents      <- frag_out$agents
  agents_grid <- frag_out$agents_grid
  
  if (switch$print_agents == 1) {
    print("FRAGMENTATION")
  }
  
  # Record immediately after fragmentation
  if ("post_fragmentation" %in% record_steps) {
    state_list[["post_fragmentation"]] <- record_state(
      step_label  = "post_fragmentation",
      step_number = steps_1 + 1
    )
  }
  

  # Post-fragmentation ------------------------------------------------------

  for (j in 1:steps_2) {
    
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
