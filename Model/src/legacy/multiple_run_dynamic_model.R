################# Multiple Runs Dynamic Model Function #######################################

# The function is another version of the dynamic model, designed to run the
# model multiple times with varying parameters and record the results into 1
# output file

multiple_run_dynamic_model <- function(steps = mod_par$steps, file_name) {

  # for displaying run time of the function
  start.time <- Sys.time()

  # create a vector with al species id's
  species_sequence <- seq(1:mod_par$n_species)

  # create a df for all constant parameter/results of the simulation
  static_output <- data.frame(
    sim_id = numeric(),
    grid_size = numeric(),
    habitat_percent = numeric(),
    spatial_ac = numeric(),
    frag_factor = numeric(),
    n_pop = numeric(),
    species = numeric(),
    niche_breadth = numeric(),
    sim_steps = numeric(),
    patches = numeric()
  )

  # create a data frame for species abundance and diversity results
  output_all <- data.frame(
    sim_id = numeric(),
    step = numeric(),
    individuals = numeric(),
    present_species = numeric(),
    shannon = numeric(),
    patch_individuals_mean = numeric(),
    patch_species_mean = numeric(),
    patch_shannon_mean = numeric(),
    scenario = character(),
    habitat = numeric(),
    fragmentation = numeric()
  )

  # Go through all parameter variation in var_par file
  for (k in 1:nrow(var_par)) {

    # Function that runs the following 3 functions: 'generate_grid', 'generate_agents', 'distribute_agents'
    model_start <- initialize(var_par$frag[k], var_par$hab[k])

    # extract simulation space, agents grid and agents list from the results list
    grid <- model_start$grid
    agents_grid <- model_start$agents_grid
    agents <- model_start$agents

    # clumping raster to get the number of patches later
    clumped <- clump(grid, direction = 4, gaps = F)
    # creating a frequency table from the patches raster
    clump_freq <- na.omit(freq(clumped))

    # create a temp df for all constant parameter/results of the simulation
    temp_static_output <- data.frame(
      sim_id = k,
      grid_size = mod_par$grid_size,
      habitat_percent = var_par$hab[k],
      spatial_ac = mod_par$spatial_ac,
      frag_factor = var_par$frag[k],
      n_pop = mod_par$n_pop,
      species = mod_par$n_species,
      niche_breadth = mod_par$niche_breadth,
      sim_steps = mod_par$steps,
      patches = clumped@data@max
    )

    # create a temp df for species abundance and diversity results
    temp_output_all <- data.frame(
      sim_id = k,
      step = rep(NA, steps),
      individuals = NA,
      present_species = NA,
      shannon = NA,
      patch_individuals_mean = NA,
      patch_species_mean = NA,
      patch_shannon_mean = NA,
      scenario = NA,
      habitat = NA,
      fragmentation = NA
    )
    # looping through the simulation time steps
    for (i in 1:steps) {

      # function that returns a vector with species presence data for the landscape
      species_count_vec <- sapply(
        species_sequence, function(sp) {
          sum(agents$species_id == sp)
        }
      )

      # creating a vector with abundance value for each patch
      abundance_patch <- sapply(clump_freq[, 1], function(id) {
        sum(agents$patch_id == id)
      })

      # creating a vector with number of species for each patch
      n_spec_patch <- sapply(clump_freq[, 1], function(id) {
        length(unique(agents$species_id[agents$patch_id == id]))
      })

      # calculating average diversity across all patches
      # a function to create a diversity vectors with diversity values of all patches
      diversity_vector <- sapply(clump_freq[, 1], function(cl) {
        single_patch <- sapply(species_sequence, function(sp) {
          sum(agents$species_id[agents$patch_id == cl] == sp)
        })

        vegan::diversity(single_patch)
      })

      # recording data in the results df
      temp_output_all$step[i] <- i
      temp_output_all$individuals[i] <- nrow(agents)
      temp_output_all$present_species[i] <- length(unique(agents$species_id))
      temp_output_all$shannon[i] <- vegan::diversity(species_count_vec, index = "shannon")
      temp_output_all$patch_individuals_mean[i] <- mean(abundance_patch)
      temp_output_all$patch_species_mean[i] <- mean(n_spec_patch)
      temp_output_all$patch_shannon_mean[i] <- mean(diversity_vector)
      temp_output_all$scenario[i] <- paste("hab = ", var_par$hab[k], " | frag = ", var_par$frag[k], sep = "")
      temp_output_all$habitat[i] <- var_par$hab[k]
      temp_output_all$fragmentation[i] <- var_par$frag[k]

      if (nrow(agents) > 0) {
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
      } else {
        print("ALL DEAD :(")
      }
    }
    static_output <- rbind(static_output, temp_static_output)
    output_all <- rbind(output_all, temp_output_all)
    print(paste(k, " simulation is completed", sep = ""))
  }

  # Save all data frames as a csv to the Outputs folder
  write.csv(output_all, paste("Outputs/", file_name, "_output_general.csv", sep = ""), row.names = F)
  write.csv(static_output, paste("Outputs/", file_name, "_simulation_info.csv", sep = ""), row.names = F)

  list_of_outputs <- list(output_all = output_all, static_output = static_output)

  end.time <- Sys.time()
  time.taken <- round(end.time - start.time, 2)
  print(time.taken)

  return(list_of_outputs)
}
