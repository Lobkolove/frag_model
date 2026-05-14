########################## Agents Generating Function ##########################

# The function returns a data.frame object listing all agents, their species
# ID number, and their environmental value which represents the optimal environmental niche.

generate_agents <- function(n_species, pop, grid, seed = NULL) {
  # mean_grid <- mean(grid@data@values[!is.na(grid@data@values)])
  # sd_grid <- sd(grid@data@values[!is.na(grid@data@values)] + 0.05)

  # Commented code generates a random environmental value while active code generates values with equal intervals

  # n_values <- runif(n = species, min = 0, max = 1) # Generate a vector of possible N values according to number of species
  # n_values <- rnorm(n_species, mean_grid, sd_grid)

  # create data frame and populating it
  withr::with_seed(seed, {
    agents_df <- data.table(
      ID = 1:pop,
      species_id = sample(1:n_species, pop, replace = T),
      x_loc = 0,
      y_loc = 0
    )
  })

  # for optimization purposes, trying to run th model with 'agents' as a matrix instead of df
  # agents_df1 <- as.matrix(agents_df)
  
  return(agents_df)
}
