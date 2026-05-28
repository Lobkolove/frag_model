######################## Birth function ########################################

# The function represent 1 of the 2 core processes in the dynamic model simulation.
# It generates new agents according to a set of rules.
# The function returns the updated agents list and agent raster layer.

birth <- function(agents, grid, nb, disp, d_dis) {

  # looping through all agents
  for (i in 1:nrow(agents)) {

    # generate random numbers between 0-1. (Should the numbers be instead pulled from a distribution?)
    rand <- runif(3, 0, 1)

    # extract N value,niche breadth, and location of current agent

    u <- species_par$n_value[species_par$species_id == agents$species_id[i]]

    # Extract species specific parameters for species of current agent
    cur_spec_par <- species_par[species_par$species_id == agents$species_id[i], ]

    # check to see if species specific niche breadth OR the function argument which comes at the moment from var_par

    if (switch$species_specific_par == 1) {
      nb <- cur_spec_par$niche_breadth
    } else if (switch$species_specific_par == 0) {
      nb <- nb
    } else {
      print("check switches please")
    }

    cur_loc <- c(agents$x_loc[i], agents$y_loc[i])
    # checking the birth rate parameter
    if (
      rand[1] < cur_spec_par$birth_rate
    ) {
      #get dispersal rate
      if (switch$species_specific_par == 1) {
        dr <- cur_spec_par$dispersal_rate
      } else if (switch$species_specific_par == 0) {
        dr <- disp
      } else {
        print("check switches please")
      }
      # checking if short or long dispersal
      if (rand[2] < dr) {
        ### short dispersal###

        # generate new location for birth
        # dispersal_kernel <- sample(c(-4, -3, -2, -1, 1, 2, 3, 4), 2)
        # new_loc_short <- c(cur_loc[1] + dispersal_kernel[1], cur_loc[2] + dispersal_kernel[2])
        new_loc_short <- toroidal_disperse(
          cur_loc = cur_loc,
          d_sd = d_dis,
          d_mean = d_dis,
          grid_size = mod_par$grid_size
        )
        inter_cell_subset <- collapse::fsubset(
          agents,
          x_loc == new_loc_short[1] & y_loc == new_loc_short[2]
        )
        intra_cell <- collapse::fnrow(collapse::fsubset(
          inter_cell_subset,
          species_id == agents$species_id[i]
        ))
        inter_cell <- collapse::fnrow(inter_cell_subset)

        if (
          # The first checks are obsolete since toroidal disperse never leaves the grid.
          # new_loc_short[1] <= mod_par$grid_size &&
          #   new_loc_short[1] >= 1 && # checking if new location is within the grid
          #   new_loc_short[2] <= mod_par$grid_size &&
          #   new_loc_short[2] >= 1 && # checking if new location is within the grid
          !is.na(grid[new_loc_short[1], new_loc_short[2]]) && # checking if the grid cell is non-matrix (habitat)
          inter_cell < mod_par$k_inter && # checking if the cell did not reach its inter-specific carrying capacity
          intra_cell < mod_par$k_intra # checking if the cell did not reach its intra-specific carrying capacity
        ) {
          # extract environmental value at the new location and calculate survival pro. according to Gravel 2006
          e <- grid[new_loc_short[1], new_loc_short[2]]
          survival_prob <- exp((-(e - u)^2) / (2 * nb^2))

          # checking the survival probability
          if (survival_prob > rand[3]) {

            # update agents list
            agents <- rbindlist(list(
              agents,
              list(
                i,
                agents$species_id[i],
                new_loc_short[1],
                new_loc_short[2]
              )
            ))
          }
        }
      } else {
        ### long dispersal###
        # generate location and calculate survival prob.
        new_loc_long <- round(runif(2, 1, mod_par$grid_size))
        e <- grid[new_loc_long[1], new_loc_long[2]]
        survival_prob <- exp((-(e - u)^2) / (2 * nb^2))

        inter_cell_subset <- collapse::fsubset(
          agents,
          x_loc == new_loc_long[1] & y_loc == new_loc_long[2]
        )
        inter_cell <- collapse::fnrow(inter_cell_subset)
        intra_cell <- collapse::fnrow(collapse::fsubset(
          inter_cell_subset,
          species_id == agents$species_id[i]
        ))

        if (
          !is.na(grid[new_loc_long[1], new_loc_long[2]]) && # checking if the grid cell is non-matrix (habitat)
            inter_cell < mod_par$k_inter && # checking if the cell did not reach its inter-specific carrying capacity
            intra_cell < mod_par$k_intra # checking if the cell did not reach its intra-specific carrying capacity
        ) {
          # checking survival prob.
          if (survival_prob > rand[3]) {

            # update agents list
            agents <- rbindlist(list(
              agents,
              list(
                i,
                agents$species_id[i],
                new_loc_long[1],
                new_loc_long[2]
              )
            ))
          }
        }
      }
    }
  }

  return(agents)
}
