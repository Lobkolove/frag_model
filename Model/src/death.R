######################## Death function ########################################

# The function represent 1 of the 2 core processes in the dynamic model simulation.
# It eliminates  agents according to a set of rules.
# The function returns the updated agents list and grid.
# Commented code (line 23, 24, 27) can consider survival prob in death process.
# The problem is that most agents already have high survival prob. and don't die.

death <- function(
  agents, 
  grid, 
  species_par,
  edge_effects = FALSE,
  edge_cells = NULL,
  species_specific = FALSE,
  death_rate = NULL,
  edge_factor = NULL
) {

  if (!species_specific && (is.null(death_rate) || is.null(edge_factor))) {
    stop("If species specific parameters are disabled, `death_rate` and `edge_factor` must be provided.")
  }

  # empty vector to record agent numbers to be deleted
  delete_agents <- vector(length = nrow(agents))
  
  # Looping through all agents
  for (i in 1:nrow(agents)) {

    # generate random values between 0-1
    rand1 <- runif(1, 0, 1)

    # Extracting agent's location, death rate and edge_effect
    cur_loc <- c(agents$x_loc[i], agents$y_loc[i])
    cur_spec_par <- species_par[species_par$species_id == agents$species_id[i], ]

    if (species_specific) {
      death_rate <- cur_spec_par$death_rate
      edge_factor <- cur_spec_par$edge_effect
    }
    
    # Updating death rate if switch is on and agent is at edge
    if (edge_effects) {
      if (is.null(edge_cells)) edge_cells <- 
      if (!is.null(edge_cells) && !is.na(edge_cells[cur_loc[1], cur_loc[2]])) {
        death_rate <- death_rate * edge_factor
      }
    }
    # checking if killing conditions are met according to death rate
    if (rand1 < death_rate) {
      # add agents to deletion list
      delete_agents[i] <- TRUE
    }
  }

  # delete agents from list
  agents <- agents[!delete_agents, , drop = FALSE]

  return(agents)
}