######################## Birth function ########################################

# The function represent 1 of the 3 core processes in the dynamic model simulation.
# It generates new agents according to a set of rules.
# The function returns the updated agents list.

birth <- function(
  agents,
  grid,
  k_inter,
  k_intra,
  species_par,
  species_specific = FALSE,
  nb = NULL,
  d_type = c("short_long", "random"),
  kernel = c("exponential", "log-normal"),
  disp = NULL,
  d_dis = NULL,
  habitat_cells = NULL
) {

  d_type <- match.arg(d_type)
  kernel <- match.arg(kernel)

  if (d_type == "short_long" && (is.null(d_dis))) {
    stop("For short/long dispersal, `d_dis` must be provided.")
  } 
  
  if (!species_specific) {
    if (d_type == "short_long" && is.null(disp)) stop("For short/long dispersal, either species specific parameters or `disp` are required.")
    if (is.null(nb)) stop("When species specific parameters are disabled, `nb` is required.")
  }

  if (nrow(grid) != ncol(grid)) stop("Non-quadratic grids are currently not supported.")
  grid_size <- nrow(grid)

  # looping through all agents
  for (i in 1:nrow(agents)) {

    # generate random numbers between 0-1. (Should the numbers be instead pulled from a distribution?)
    rand <- runif(3, 0, 1)
    
    # Extract species specific parameters for species of current agent
    cur_spec_id <- agents$species_id[i]
    cur_spec_par <- species_par[species_par$species_id == cur_spec_id, ]
    
    # extract N value and location of current agent

    # n value
    u <- cur_spec_par$n_value
    # parent location
    cur_loc <- c(agents$x_loc[i], agents$y_loc[i])
    
    # if species specific parameters are enabled, niche breadth and dispersal rate are taken from species_par
    # instead of 
    if (species_specific) {
      nb <- cur_spec_par$niche_breadth
      disp <- cur_spec_par$dispersal_rate
    }

    # checking the birth rate parameter
    if (rand[1] < cur_spec_par$birth_rate) {

      # checking dispersal type (random - nonrandom)
      if (d_type == "random") {

        new_loc_hab <- random_disperse(grid, habitat_cells = habitat_cells)

        # check if the new location is a habitat cell (skip if not)
        if (is.na(grid[new_loc_hab[1], new_loc_hab[2]])) {
          stop("Random dispersal sent to a non-habitat cell. Please check the `habitat_cells` argument.")
        }
        
        # extract environmental value at the new location and calculate survival pro. according to Gravel 2006
        e <- grid[new_loc_hab[1], new_loc_hab[2]]
        survival_prob <- exp((-(e - u)^2) / (2 * nb^2))

        # checking the survival probability
        if (survival_prob > rand[3]) {
          inter_cell_subset <- collapse::fsubset(
            agents,
            x_loc == new_loc_hab[1] & y_loc == new_loc_hab[2]
          )
          intra_cell <- collapse::fnrow(collapse::fsubset(
            inter_cell_subset,
            species_id == cur_spec_id
          ))
          inter_cell <- collapse::fnrow(inter_cell_subset)

          if (
            inter_cell < k_inter && # checking if the cell did not reach its inter-specific carrying capacity
              intra_cell < k_intra # checking if the cell did not reach its intra-specific carrying capacity
          ) {
            # update agents list
            new_row <- data.table(
              ID = i,
              species_id = cur_spec_id,
              x_loc = new_loc_hab[1],
              y_loc = new_loc_hab[2]
            )
            agents <- rbind(agents, new_row, use.names = TRUE)
          }
        }
      } else if (d_type == "short_long") {
        
        # checking if short or long dispersal
        if (rand[2] < disp) {
          ### short dispersal###
  
          # Get new location
          new_loc_short <- toroidal_disperse(
            cur_loc = cur_loc,
            d_sd = d_dis,
            d_mean = d_dis,
            grid_size = grid_size,
            kernel = kernel
          )
  
          if (!is.na(grid[new_loc_short[1], new_loc_short[2]])) {
            # checking if the grid cell is non-matrix (habitat)
            # extract environmental value at the new location and calculate survival pro. according to Gravel 2006
            e <- grid[new_loc_short[1], new_loc_short[2]]
            survival_prob <- exp((-(e - u)^2) / (2 * nb^2))
  
            # checking the survival probability
            if (survival_prob > rand[3]) {
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
                inter_cell < k_inter && # checking if the cell did not reach its inter-specific carrying capacity
                  intra_cell < k_intra # checking if the cell did not reach its intra-specific carrying capacity
              ) {
                # update agents list
                new_row <- data.table(
                  ID = i,
                  species_id = cur_spec_id,
                  x_loc = new_loc_short[1],
                  y_loc = new_loc_short[2]
                )
                agents <- rbind(agents, new_row, use.names = TRUE)
              }
            }
          }
        } else {
          ### long dispersal###
          # generate location and calculate survival prob.
          new_loc_long <- random_disperse(
            grid,
            force_habitat = FALSE
          )

          # check if the new location is a habitat cell (skip if not)
          if (is.na(grid[new_loc_long[1], new_loc_long[2]])) next 
          
          # extract environmental value at the new location and calculate survival pro. according to Gravel 2006
          e <- grid[new_loc_long[1], new_loc_long[2]]
          survival_prob <- exp((-(e - u)^2) / (2 * nb^2))

          # check survival probability (skip if didn't survive)
          if (survival_prob <= rand[3]) next
  
          inter_cell_subset <- collapse::fsubset(
            agents,
            x_loc == new_loc_long[1] & y_loc == new_loc_long[2]
          )
          inter_cell <- collapse::fnrow(inter_cell_subset)
          intra_cell <- collapse::fnrow(collapse::fsubset(
            inter_cell_subset,
            species_id == cur_spec_id
          ))
  
          if (
              inter_cell < k_inter && # checking if the cell did not reach its inter-specific carrying capacity
              intra_cell < k_intra # checking if the cell did not reach its intra-specific carrying capacity
          ) {
              # update agents list
              new_row <- data.table(
                ID = i,
                species_id = cur_spec_id,
                x_loc = new_loc_long[1],
                y_loc = new_loc_long[2]
              )
              agents <- rbind(agents, new_row, use.names = TRUE)
            }
        }
      }
    }
  }

  return(agents)
}
