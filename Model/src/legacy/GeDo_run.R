################# GeDo run Dynamic Model Function ##############################

# The function is a third iteration of the dynamic model, designed to separate
# geometric and demographic fragmentation effects by following these steps:
# 1. run the model with 100% habitat
# 2. apply the cookie-cutter approach by fragmenting the landscape (the effects at this time point are geometric)
# 3. keep running model and observe the demographic effects

# Functions to initialize empty result files

source("Model/src/initialize_result_files.R")

GeDo_run <- function(mod_par, 
                     var_par, 
                     switch,
                     file_name = "test",
                     task_id = "no_cluster",
                     sim_id) {
  # create a vector with all species id's
  species_sequence <- 1:mod_par$n_species
  steps_1 <- mod_par$steps_pre_frag
  steps_2 <- mod_par$steps_post_frag


  # Initialize output -------------------------------------------------------

  # create a data frame for species abundance and diversity results in the whole landscape
  output_all <- intialize_output_all()
  # create a data frame for species abundance and diversity results in the whole landscape
  output_sample <- initialize_output_sample(species_seq = species_sequence)


  # Run model ---------------------------------------------------------------

  # Go through all parameter variation in var_par data frame
  start.time.sim <- Sys.time()

  # either stochastic or not-stochastic spin-up phase. Can be switched in parameter file
  if (switch$random == 1) {
    set.seed(seed)
  }

  # initialize the model with 100% habitat
  model_start <- initialize(
    frag = 0,
    hab = 1,
    ac = var_par$ac,
    nb = var_par$nb
  )

  # extract simulation space, agents grid and agents list from the results list
  grid <- model_start$grid
  agents_grid <- model_start$agents_grid
  agents <- model_start$agents

  # create a temp df for species abundance and diversity results

  # looping through the simulation time steps
  for (i in 1:steps_1) {
    start.time <- Sys.time()

    # function that returns a vector with species presence data for the landscape
    species_count_vec <- sapply(
      species_sequence, function(sp) {
        sum(agents$species_id == sp)
      }
    )

    # recording data in the results df
    output_all <- tibble::add_row(output_all,
      sim_id = sim_id,
      step = i,
      individuals = nrow(agents),
      present_species = length(unique(agents$species_id)),
      shannon = vegan::diversity(species_count_vec, index = "shannon"),
      patch_individuals_mean = nrow(agents),
      patch_species_mean = length(unique(agents$species_id)),
      patch_shannon_mean = vegan::diversity(species_count_vec, index = "shannon"),
      sample_individuals_mean = 0,
      sample_species_mean = 0,
      sample_shannon_mean = 0,
      habitat = 1,
      fragmentation = 0,
      n_patches = 1,
      n_samples = 0
    )

    if (nrow(agents) > 0) {
      # run birth function and update agents list and grid
      step1 <- birth(
        agents = agents,
        agents_grid = agents_grid,
        grid = grid,
        NB = var_par$nb,
        disp = var_par$disp,
        d_dis = var_par$disp_dist
      )
      agents <- step1$agents
      agents_grid <- step1$agents_grid

      # run death function and update agents list and grid
      step2 <- death(
        agents = agents,
        agents_grid = agents_grid,
        grid = grid,
        edge_fac = var_par$edge
      )
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

    end.time <- Sys.time()
    time.taken <- round(end.time - start.time, 2)
    if (switch$print_agents == 1) {
      print(paste("step ", i, " took ", time.taken, " with ", nrow(agents), " agents", sep = ""))
    }
    if (switch$animation_export == 1) { # export png to animation folder in case switch is turned on
      png(paste("Outputs/animation/", sprintf("%04d", i), "_animation_plot.png", sep = ""), width = 600, height = 640, pointsize = 10)
      plot(t(flip(grid, 1)), col = gray.colors(21), legend = FALSE, axes = FALSE) # plotting the simulation space
      plot(t(flip(agents_grid, 1)), breaks = breaks_agents, col = rainbow(mod_par$n_species), add = T, axes = F) # adding the agents raster
      plot(t(flip(grid, 1)), col = gray.colors(21), breaks = breaks_space, legend.only = T, horizontal = T, ) # adding legend of the simulation space
      text(x = 0, y = mod_par$grid_size + mod_par$grid_size / 10, adj = c(0, 1), labels = paste("Time Step = ", i, sep = ""), cex = 2, xpd = T)

      dev.off()
    }
  }
  # randomize again for fragmentation
  if (switch$random_post_frag == 0) {
    set.seed(NULL)
  } else {
    set.seed(seed)
  }

  # cookie-cut the landscape
  fragmanted_space <- cookie_cutting(grid, agents, agents_grid, var_par$hab, var_par$frag)

  grid <- fragmanted_space$grid
  agents <- fragmanted_space$agents
  agents_grid <- fragmanted_space$agents_grid

  if (switch$print_agents == 1) {
    print("FRAGMENTATION")
  }

  # choose habitat cells to sample
  grid_values <- getValues(grid)
  possible_cells <- which(!is.na(grid_values))

  # choose between all cells or number of samples according to sample parameter
  if (switch$sample_all == 1) {
    samples <- possible_cells
  } else if (switch$sample_all == 0) {
    samples <- sample(possible_cells, mod_par$n_samples)
  } else {
    print("check switches please")
  }
  # clumping raster to get the number of patches later
  clumped <- clump(grid, directions = 4, gaps = F)
  # creating a frequency table from the patches raster
  clump_freq <- na.omit(freq(clumped))

  for (j in 1:steps_2) {
    start.time <- Sys.time()
    # an if statement to record data only at last time step for optimization - delete when all data needed

    # if (j == steps_2) { # a loop to increase simulation speed in case results are needed only for the last time step

    # function that returns a vector with species presence data for the landscape
    species_count_vec <- sapply(
      species_sequence, function(sp) {
        sum(agents$species_id == sp)
      }
    )

    # creating a vector with abundance value for each patch
    abundance_patch <- sapply(clump_freq[, 1], function(id) {
      sum(agents$patch_id == id, na.rm = T)
    })

    # creating a vector with number of species for each patch
    n_spec_patch <- sapply(clump_freq[, 1], function(id) {
      length(unique(agents$species_id[agents$patch_id == id]))
    })

    # calculating average diversity across all patches
    # a function to create a diversity vectors with diversity values of all patches
    diversity_vector <- sapply(clump_freq[, 1], function(cl) {
      single_patch <- sapply(species_sequence, function(sp) {
        sum(agents$species_id[agents$patch_id == cl] == sp, na.rm = T)
      })

      vegan::diversity(single_patch)
    })

    # vector of abundance data for all samples

    abundance_sample <- sapply(samples, function(sm) {
      xyloc <- rowColFromCell(grid, sm)
      sum(agents$x_loc == xyloc[1] & agents$y_loc == xyloc[2])
    })

    # vector for richness data for all samples

    richness_sample <- sapply(samples, function(sm) {
      xyloc <- rowColFromCell(grid, sm)
      length(unique(agents$species_id[agents$x_loc == xyloc[1] & agents$y_loc == xyloc[2]]))
    })

    # vector of shannon diversity for each sample

    diversity_vector_sample <- sapply(samples, function(sm) {
      single_sample <- sapply(species_sequence, function(sp) {
        xyloc <- rowColFromCell(grid, sm)
        sum(agents$species_id[agents$x_loc == xyloc[1] & agents$y_loc == xyloc[2]] == sp, na.rm = T)
      })

      vegan::diversity(single_sample)
    })

    # recording data in the results df
    output_all <- tibble::add_row(output_all,
                                  sim_id = sim_id,
                                  step = j + steps_1,
                                  individuals = nrow(agents),
                                  present_species = length(unique(agents$species_id)),
                                  shannon = vegan::diversity(species_count_vec, index = "shannon"),
                                  patch_individuals_mean = mean(abundance_patch),
                                  patch_species_mean = mean(n_spec_patch),
                                  patch_shannon_mean = mean(diversity_vector),
                                  sample_individuals_mean = mean(abundance_sample),
                                  sample_species_mean = mean(richness_sample),
                                  sample_shannon_mean = mean(diversity_vector_sample),
                                  habitat = var_par$hab,
                                  fragmentation = var_par$frag,
                                  n_patches = clumped@data@max,
                                  n_samples = length(samples)
    )
    
    # recording data in the sample output data frame at first time step after fragmentation and at the last time step
    if (j == 1 | j == steps_2) {
      # go through all samples 1 by 1
      for (l in 1:length(samples)) {
        xyloc <- rowColFromCell(grid, samples[l]) # get sample coordinates
        sample_patch_id <- clumped[samples[l]] # get the sample patch id

        # create a presence vector for sample l and for all species
        species_presence_vec <- sapply(
          species_sequence, function(sp) {
            sum(agents$species_id == sp & agents$x_loc == xyloc[1] & agents$y_loc == xyloc[2])
          }
        )
        # add all data that will go into "output_sample" in 1 row as a vector. Reason for vector instead of df row is the changing number of species between simulations
        new_row <- c(sim_id, 
                     clumped[xyloc], 
                     l, 
                     j + steps_1, 
                     xyloc[1], 
                     xyloc[2], 
                     clump_freq[sample_patch_id, 2], 
                     var_par$frag, 
                     var_par$hab, 
                     species_presence_vec)

        output_sample <- rbind(output_sample, new_row)
      }
    }
    if (nrow(agents) > 0) {
      # run birth function and update agents list and grid
      step1 <- birth(
        agents = agents,
        agents_grid = agents_grid,
        grid = grid,
        NB = var_par$nb,
        disp = var_par$disp,
        d_dis = var_par$disp_dist
      )
      agents <- step1$agents
      agents_grid <- step1$agents_grid

      # run death function and update agents list and grid
      step2 <- death(
        agents = agents,
        agents_grid = agents_grid,
        grid = grid,
        edge_fac = var_par$edge
      )
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
    n_agents <- nrow(agents)
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time, 2)
    if (switch$print_agents == 1) {
      print(paste("step ", j + steps_1, " took ", time.taken, " with ", n_agents, " agents", sep = ""))
    }

    if (switch$animation_export == 1) { # export png to animation folder in case switch is turned on
      png(paste("Outputs/animation/", sprintf("%04d", j + i), "_animation_plot.png", sep = ""), width = 600, height = 640, pointsize = 10)
      plot(t(flip(grid, 1)), col = gray.colors(21), legend = FALSE, axes = FALSE) # plotting the simulation space
      plot(t(flip(agents_grid, 1)), breaks = breaks_agents, col = rainbow(mod_par$n_species), add = T, axes = F) # adding the agents raster
      plot(t(flip(grid, 1)), col = gray.colors(21), breaks = breaks_space, legend.only = T, horizontal = T, ) # adding legend of the simulation space
      text(x = 0, y = mod_par$grid_size + mod_par$grid_size / 10, adj = c(0, 1), labels = paste("Time Step = ", j + i, sep = ""), cex = 2, xpd = T)

      dev.off()
    }
  }

  print(paste(" simulation ", sim_id, " is completed", sep = ""))

  end.time.sim <- Sys.time()
  time.taken.sim <- round(end.time.sim - start.time.sim, 2)
  print(paste(" and it took ", time.taken.sim, sep = ""))

  # prepare grid raster for export

  grid_sam <- grid
  grid_sam[grid_sam] <- NA
  grid_sam[samples] <- 2

  # save raster of simulation space and of samples

  if (switch$export_raster == 1) {
    writeRaster(grid, paste("Outputs/", file_name, "_sim_", sim_id, "_sim_space_raster.grd", sep = ""))
    writeRaster(grid_sam, paste("Outputs/", file_name, "_sim_", sim_id, "_samples_raster.grd", sep = ""))
  } else if (switch$export_raster == 0) {

  } else {
    print("check switches please")
  }


  # Write output ------------------------------------------------------------

  # remove 1st row of samples df
  output_sample <- output_sample[-1, ]

  list_of_outputs <- list(output_all = output_all, output_sample = output_sample)

  return(list_of_outputs)
}
