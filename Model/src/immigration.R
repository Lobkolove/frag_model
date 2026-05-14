######################## Immigration function ##################################

# The function is an additional process in the model which allows immigration of
# individuals from the species pool back into the simulation.
# This allows "extinct" species to re-establish in the system.

immigration <- function(agents, grid) {

  for (i in 1:mod_par$n_immigrants) {
    rand <- runif(1, 0, 1) # random number for survival prob
    cur_spec <- sample(1:mod_par$n_species, 1) # sample a random species from species pool
    rand_loc <- round(runif(2, 1, mod_par$grid_size)) # choose a random location on the grid
    u <- species_par$n_value[species_par$species_id == cur_spec] # get environmental optimum of current species
    intra_cell <- nrow(agents[agents$x_loc == rand_loc[1] & agents$y_loc == rand_loc[2] & agents$species_id == agents$species_id[i], ]) # get the number of same species individuals in the random cell
    inter_cell <- nrow(agents[agents$x_loc == rand_loc[1] & agents$y_loc == rand_loc[2], ]) # get number of individuals in the cell
    nb <- species_par$niche_breadth[species_par$species_id == cur_spec]
    
    if (!is.na(grid[rand_loc[1], rand_loc[2]]) && # check the location is non-matrix
      inter_cell < mod_par$k_inter && # check the cell didn't reach inter specific capacity
      intra_cell < mod_par$k_intra) { # check the cell didn't reach intra specific capacity

      e <- grid[rand_loc[1], rand_loc[2]] # extract the n value of the cell
      survival_prob <- exp((-(e - u)^2) / (2 * nb^2)) # calculating survival prob

      if (survival_prob > rand) {

        # update agents list
        new_row <- data.table(
          ID = i,
          species_id = cur_spec,
          x_loc = rand_loc[1],
          y_loc = rand_loc[2]
        )
        agents <- rbindlist(list(agents, new_row), use.names = T)
      }
    }
  }

  return(agents)
}