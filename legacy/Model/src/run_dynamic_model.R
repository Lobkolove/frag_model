################# Dynamic model function #######################################

# The function initialize the model and then goes through birth and death as
# many times as specified by "steps". 3 output files are created and exported to
# the "Output" folder with a unique file name depending on the 2nd argument to the
# function -  "file_name". The function returns the 3 output files as a list

run_dynamic_model <- function(frag, hab, steps = mod_par$steps, file_name) {

  # for displaying run time of the function
  start.time <- Sys.time()

  # set breaks for visualization
  breaks_agents <- seq(1, mod_par$n_species)
  breaks_space <- seq(0, 1, by = 0.05)

  # Function that runs the following 3 functions: 'generate_grid', 'generate_agents', 'distribute_agents'
  model_start <- initialize(frag, hab)

  # extract simulation space, agents grid and agents list from the results list
  grid <- model_start$grid
  agents_grid <- model_start$agents_grid
  agents <- model_start$agents

  # clumping raster to get the number of patches later
  clumped <- clump(grid, direction = 4)

  # create a df for all constant parameter/results of the simulation
  static_output <- data.frame(
    grid_size = mod_par$grid_size,
    habitat_percent = mod_par$habitat_percent,
    spatial_ac = mod_par$spatial_ac,
    frag_factor = mod_par$frag_factor,
    n_pop = mod_par$n_pop,
    species = mod_par$n_species,
    niche_breadth = mod_par$niche_breadth,
    sim_steps = mod_par$steps,
    patches = clumped@data@max
  )

  # create a df for species abundance and diversity results
  output_all <- data.frame(
    step = rep(NA, steps),
    individuals = NA,
    present_species = NA,
    shannon = NA
  )

  # a table for recording species presence
  output_species <- data.frame(matrix(nrow = 0, ncol = mod_par$n_species + 1))

  # looping through the simulation time steps
  for (i in 1:steps) {

    # function that returns a vector with species presence data
    species_sequence <- seq(1:mod_par$n_species)
    species_count_vec <- sapply(
      species_sequence, function(sp) {
        sum(agents$species_id == sp)
      }
    )

    # recording data in the results df
    output_all$step[i] <- i
    output_all$individuals[i] <- nrow(agents)
    output_all$present_species[i] <- length(unique(agents$species_id))
    output_all$shannon[i] <- vegan::diversity(species_count_vec, index = "shannon")
    output_species <- rbind(output_species, c(i, species_count_vec))

    # run birth function and update agents list and grid
    step1 <- birth(agents, agents_grid, grid)
    agents <- step1$agents
    agents_grid <- step1$agents_grid

    # run death function and update agents list and grid
    step2 <- death(agents, agents_grid, grid)
    agents <- step2$agents
    agents_grid <- step2$agents_grid

    # run immigration function and update agent list and grid
    if (switch$immigration == 1) {
      step3 <- immigration(agents, agents_grid, grid)
      agents <- step3$agents
      agents_grid <- step3$agents_grid
    }

    print(paste("abundance = ", nrow(agents), sep = ""))
    print(paste("sp richness = ", length(unique(agents$species_id)), sep = ""))

    # save a plot of the time step for later animation (can be switched off)
    if (switch$animation_export == 1) {
      png(paste("Outputs/animation/", sprintf("%04d", i), "_animation_plot.png", sep = ""), width = 600, height = 640, pointsize = 10)
      plot(t(flip(grid, 1)), col = gray.colors(21), legend = F, axes = F) # plotting the simulation space
      plot(t(flip(agents_grid, 1)), breaks = breaks_agents, col = rainbow(mod_par$n_species), add = T, axes = F) # adding the agents raster
      plot(t(flip(grid, 1)), col = gray.colors(21), breaks = breaks_space, legend.only = T, horizontal = T, ) # adding legend of the simulation space
      text(x = 0, y = mod_par$grid_size + mod_par$grid_size / 12, adj = c(0, 1), labels = paste("Time Step = ", i, sep = ""), cex = 2, xpd = T)

      dev.off()
    }
  }

  # renaming species df (for some reason colnames get deleted if done earlier)
  colnames(output_species) <- c("step", paste("specie_", 1:mod_par$n_species, sep = ""))

  # Save all data frames as a csv to the Outputs folder
  write.csv(output_all, paste("Outputs/", file_name, "_output_general.csv", sep = ""), row.names = F)
  write.csv(output_species, paste("Outputs/", file_name, "_output_species.csv", sep = ""), row.names = F)
  write.csv(static_output, paste("Outputs/", file_name, "_simulation_info.csv", sep = ""), row.names = F)

  list_of_outputs <- list(output_all = output_all, static_output = static_output, output_species = output_species, grid = grid, agents = agents, agents_grid = agents_grid)

  end.time <- Sys.time()
  time.taken <- round(end.time - start.time, 2)
  print(time.taken)

  return(list_of_outputs)
}