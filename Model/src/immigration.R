######################## Immigration function ##################################

# The function is an additional process in the model which allows immigration of
# individuals from the species pool back into the simulation.
# This allows "extinct" species to re-establish in the system.

immigration <- function(
  agents,
  grid,
  n_immigrants,
  n_species,
  k_inter,
  k_intra,
  species_par,
  species_specific = FALSE,
  nb = NULL
) {

  if (!species_specific && is.null(nb)) {
    stop("If species specific parameters are disabled, `nb` must be provided.")
  } else if (species_specific && !is.null(nb)) {
    warning("Species specific parameters enabled: ignoring `nb`.")
  }

  if (nrow(grid) != ncol(grid)) {
    stop("Non-quadratic grids are currently not supported.")
  }
  grid_size <- nrow(grid)

  for (i in 1:n_immigrants) {
    rand <- runif(1, 0, 1) # random number for survival prob
    cur_spec <- sample(1:n_species, 1) # sample a random species from species pool
    rand_loc <- round(runif(2, 1, grid_size)) # choose a random location on the grid
    if (is.na(grid[rand_loc[1], rand_loc[2]])) {
      next # skip if the location is a matrix
    }

    u <- species_par$n_value[species_par$species_id == cur_spec] # get environmental optimum of current species
    intra_cell <- nrow(agents[
      agents$x_loc == rand_loc[1] &
        agents$y_loc == rand_loc[2] &
        agents$species_id == cur_spec,
    ]) # get number of same species individuals in the chosen cell
    inter_cell <- nrow(agents[
      agents$x_loc == rand_loc[1] & agents$y_loc == rand_loc[2],
    ]) # get total number of individuals in the cell
    
    if (
      inter_cell < k_inter && # check the cell didn't reach inter specific capacity
        intra_cell < k_intra # check the cell didn't reach intra specific capacity
    ) {
      
      if (species_specific) {
        nb <- species_par$niche_breadth[species_par$species_id == cur_spec]
      } 
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
        agents <- rbind(agents, new_row, use.names = TRUE)
      }
    }
  }

  return(agents)
}