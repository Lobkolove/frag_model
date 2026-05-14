####################### Agent Distribution Function ############################

# The function returns a raster layer where cells values are equal to the species ID number
# indicating where agents are distributed in space.

distribute_agent <- function(gr_size, 
                             agents, 
                             space, 
                             nb,
                             seed = NULL,
                             random_distribution = FALSE,
                             ...) {

  # extract the extent of the simulation space raster, convert to matrix
  extent <- extent(space)
  space_mx <- as.matrix(space)
  grid_values <- getValues(space)

  samp <- function(x, ...) x[sample.int(length(x), ...)] # redefine sample to work in case of vector length 1

  # populate matrix if random location is a habitat and survival probability is higher than rand1
  for (k in 1:length(agents$ID)) {

    # extract possible locations for the species 
    #nb <- species_par$niche_breadth[species_par$species_id == agents$species_id[k]] #commented out since it overrides the nb value from initialize call
    n_value <- species_par$n_value[species_par$species_id == agents$species_id[k]]
    range <- c(n_value - nb, n_value + nb)
    possible_cells <- which(grid_values > range[1] & grid_values < range[2])

    # Added to allow for random distribution of agents across the landscape, regardless of their niche values. 
    # This is useful for testing the effects of fragmentation without the confounding effect of niche-based distribution.
    if (random_distribution) {
      possible_cells <- which(!is.na(grid_values))
    }

    if (length(possible_cells) >= 1) {
      location <- samp(possible_cells, 1)
      xyloc <- rowColFromCell(space, location) # get the coordinates of the new location

      if (
        collapse::fnrow(
          collapse::fsubset(agents, x_loc == xyloc[1] & # checking if the cell did not reach its inter specific cell capacity
                    y_loc == xyloc[2])) < mod_par$k_inter &&
        collapse::fnrow(
          collapse::fsubset(agents, x_loc == xyloc[1] & # checking if the cell did not reach its intra specific cell capacity
          y_loc == xyloc[2] &
          species_id == agents$species_id[k])) < mod_par$k_intra
      ) {
        agents$x_loc[k] <- xyloc[1]
        agents$y_loc[k] <- xyloc[2]
      }
    }
  }
  
  # updating agents list
  agents <- agents[agents$x_loc != 0, ]

  return(agents)
  }
