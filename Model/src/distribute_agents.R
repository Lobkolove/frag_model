####################### Agent Distribution Function ############################

# The function returns a raster layer where cells values are equal to the species ID number
# indicating where agents are distributed in space.

distribute_agent <- function(gr_size, 
                             agents_list = agents, 
                             space, 
                             nb,
                             seed = NULL,
                             random_distribution = FALSE,
                             ...) {

  # extract the extent of the simulation space raster, convert to matrix and create blank matrix
  extent <- extent(space)
  space_mx <- as.matrix(space)
  agents_mx <- matrix(nrow = gr_size, ncol = gr_size)
  clumped <- clump(space, directions = 4)
  samp <- function(x, ...) x[sample.int(length(x), ...)] # redefine sample to work in case of vector length 1
  grid_values <- getValues(space)


  # populate matrix if random location is a habitat and  survival probability is higher than rand1
  for (k in 1:length(agents_list$ID)) {

    # extract possible locations for the species 
    #nb <- species_par$niche_breadth[species_par$species_id == agents_list$species_id[k]] #commented out since it overrides the nb value from initialize call
    n_value <- species_par$n_value[species_par$species_id == agents_list$species_id[k]]
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
          collapse::fsubset(agents_list, x_loc == xyloc[1] & # checking if the cell did not reach its inter specific cell capacity
                    y_loc == xyloc[2])) < mod_par$k_inter &&
        collapse::fnrow(
          collapse::fsubset(agents_list, x_loc == xyloc[1] & # checking if the cell did not reach its intra specific cell capacity
          y_loc == xyloc[2] &
          species_id == agents_list$species_id[k])) < mod_par$k_intra
      ) {
        agents_mx[xyloc[1], xyloc[2]] <- agents_list$species_id[k]
        agents_list$x_loc[k] <- xyloc[1]
        agents_list$y_loc[k] <- xyloc[2]
        agents_list$patch_id[k] <- clumped[xyloc[1], xyloc[2]]
      }
    }
  }
  
  # updating agents list
  agents_list <- agents_list[agents_list$x_loc != 0, ]

  # rasterize the agents location matrix
  agents_raster <- raster(agents_mx)

  # setting the same extent of the simulation space raster
  agents_raster <- setExtent(agents_raster, extent)

  return_list <- list(agents_raster = agents_raster, agents_list = agents_list)
  return(return_list)
  }
