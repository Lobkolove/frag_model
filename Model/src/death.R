######################## Death function ########################################

# The function represent 1 of the 2 core processes in the dynamic model simulation.
# It eliminates  agents according to a set of rules.
# The function returns the updated agents list and grid.
# Commented code (line 23, 24, 27) can consider survival prob in death process.
# The problem is that most agents already have high survival prob. and don't die.

death <- function(agents, grid, edge_fac) {

  # empty vector to record agent numbers to be deleted
  delete_agents <- vector()

  # crate a grid of edges
  edge <- boundaries(grid, type = "inner", asNA = TRUE)
  
  # Looping through all agents
  for (i in 1:nrow(agents)) {

    # generate random values between 0-1 and extracting agent's N value.
    rand1 <- runif(1, 0, 1)
    rand2 <- runif(1, 0, 1)

    # Extracting agent's location and death rate
    cur_loc <- c(agents$x_loc[i], agents$y_loc[i])
    death_rate <- species_par$death_rate[species_par$species_id == agents$species_id[i]]

    #getting correct edge effect factor
    if (switch$species_specific_par == 1){
      ee <- species_par$edge_effect[species_par$species_id == agents$species_id[i]]
    } else if (switch$species_specific_par == 0){
      ee <- edge_fac
    } else {
      print("check switches please")
    }
    
    # Updating death rate if switch is on and agent is at edge
    if (switch$edge_effect == 1 && !is.na(edge[cur_loc[1], cur_loc[2]])) {
      death_rate <- death_rate * ee
    }
    # checking if killing conditions are met according to death rate
    if (rand1 < death_rate) {
      # add agents to deletion list
      delete_agents <- append(delete_agents, i)
    }
  }

  # delete agents from list
  agents <- agents[-c(delete_agents), ]

  return(agents)
}