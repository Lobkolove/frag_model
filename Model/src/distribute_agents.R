####################### Agent Distribution Function ############################

distribute_agent <- function(
  agents,
  grid,
  species_par,
  k_inter,
  k_intra,
  nb = NULL,
  species_specific = FALSE,
  ...
) {

  if (!species_specific && is.null(nb)) {
    stop("When species specific parameters are disabled, `nb` is required.")
  } else if (species_specific && !is.null(nb)) {
    warning("Species specific parameters enabled: ignoring `nb`.")
  }
  
  # extract the extent of the simulation grid raster, convert to matrix
  extent <- extent(grid)
  grid_values <- getValues(grid)

  samp <- function(x, ...) x[sample.int(length(x), ...)] # redefine sample to work in case of vector length 1

  # populate matrix if random location is a habitat and survival probability is higher than rand1
  for (i in 1:nrow(agents)) {

    cur_spec_id <- agents$species_id[i]
    cur_spec_par <- species_par[species_par$species_id == cur_spec_id, ]
    u <- cur_spec_par$n_value

    if (species_specific) {
      nb <- cur_spec_par$niche_breadth
    }

    range <- c(u - nb, u + nb)
    possible_cells <- which(grid_values > range[1] & grid_values < range[2])

    if (length(possible_cells) >= 1) {
      destination <- samp(possible_cells, 1)
      xyloc <- rowColFromCell(grid, destination) # get the coordinates of the new location

      if (
        collapse::fnrow(
          collapse::fsubset(
            agents,
            x_loc == xyloc[1] & # checking if the cell did not reach its inter specific cell capacity
              y_loc == xyloc[2]
          )
        ) <
          k_inter &&
          collapse::fnrow(
            collapse::fsubset(
              agents,
              x_loc == xyloc[1] & # checking if the cell did not reach its intra specific cell capacity
                y_loc == xyloc[2] &
                species_id == cur_spec_id
            )
          ) <
            k_intra
      ) {
        agents$x_loc[i] <- xyloc[1]
        agents$y_loc[i] <- xyloc[2]
      }
    }
  }

  # updating agents list
  agents <- agents[agents$x_loc != 0, ]

  return(agents)
}
