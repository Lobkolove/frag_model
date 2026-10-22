
########################## Multiple runs Function ##############################

# The function runs the simulation a set number of time for each combination of parameters that is specified in the 'var_par' file.
# Simulation results are recorded in 3 different tables and are returned by the function as a list.
# BUG - When running the function for the first time without a single run before the function can't seem to find the objecy 'agents'
# The problem has something to do with the environment where 'agents' live, but I could not find a solution for it so far.
# THe function does generate a new grid and a new distribution raster for each run but I can not say for sure if a new 'agents'
# data.frame is being created for each run.

multi_runs <- function(runs) {
  
  start.time <- Sys.time() # for displaying run time of the function
  
  counter <- 0 # counter is used when recording sim_id for the patches output file
  
  output_species_all <- data.frame(matrix(nrow = 0, ncol = mod_par$n_species)) # creating an empty data.frame for species output
  
  # renaming the columns of the species data.frame
  colnames(output_species_all) <- paste("species_", 1:mod_par$n_species, sep = "")
  
  output_patches_all <- data.frame(
    sim_id = numeric(), # creating an empty data.frame for the patches information
    patch_id = numeric(),
    patch_area = numeric(),
    n_species = numeric(),
    individuals = numeric()
  )
  
  output_all <- data.frame(
    sim_id = numeric(), # creating an empty data.frame for general simulation output
    grid_size = numeric(),
    habitat_percent = numeric(),
    spatial_ac = numeric(),
    frag_factor = numeric(),
    n_pop = numeric(),
    species = numeric(),
    niche_breadth = numeric(),
    total_patches = numeric(),
    individuals = numeric(),
    present_species = numeric(),
    shannon = numeric(),
    N_mean = numeric(),
    N_1st_qu = numeric(),
    N_3rd_qu = numeric()
  )
  
  for (k in 1:nrow(var_par)) { # looping through the parameter data.frame created using grid.expand()
    
    # reading the current parameter settings for fragmentation, auto-correlation, and habitat percent
    
    current_frag <- var_par[k, "frag"]
    current_ac <- var_par[k, "ac"]
    current_hab <- var_par[k, "hab"]
    
    output <- data.frame(
      sim_id = rep(0, runs), # creating a small output data.frame for general data (later adding to "output_all")
      grid_size = 0,
      habitat_percent = 0,
      spatial_ac = 0,
      frag_factor = 0,
      n_pop = 0,
      species = 0,
      niche_breadth = 0,
      total_patches = 0,
      individuals = 0,
      present_species = 0,
      shannon = 0,
      N_mean = 0,
      N_1st_qu = 0,
      N_3rd_qu = 0
    )
    
    output_species <- data.frame(matrix(nrow = runs, ncol = mod_par$n_species)) # creating a small output data.frame for species data (later adding to "output_species_all")
    
    colnames(output_species) <- paste("species_", 1:mod_par$n_species, sep = "")
    
    for (l in 1:runs) { # looping through the specified runs using a single parameter setting
      
      grid <- gen_grid(mod_par$grid_size, current_ac, current_frag, current_hab)
      agents <- gen_agent(mod_par$n_species, mod_par$n_pop)
      agents_grid <- dist_agent(mod_par$grid_size, agents, grid, mod_par$niche_breadth)
      
      counter <- counter + 1
      
      # recording data in general output data.frame
      output$grid_size[l] <- mod_par$grid_size
      output$habitat_percent[l] <- current_hab
      output$spatial_ac[l] <- current_ac
      output$frag_factor[l] <- current_frag
      output$n_pop[l] <- mod_par$n_pop
      output$species[l] <- mod_par$n_species
      output$niche_breadth[l] <- mod_par$niche_breadth
      
      clumped <- clump(grid, direction = 4) # creating a raster where cell values are patch ID
      output$total_patches[l] <- clumped@data@max # recording total number of patches
      output$individuals[l] <- sum(!is.na(agents_grid@data@values)) # recording total amount of agents
      output$present_species[l] <- length(unique(agents_grid@data@values)) - 1 # not perfect but using -1 assuming that there will always be an empty cell (NA)
      
      output$N_mean[l] <- mean(grid@data@values, na.rm = T)
      output$N_1st_qu[l] <- quantile(grid, 0.25, na.rm = T)
      output$N_3rd_qu[l] <- quantile(grid, 0.75, na.rm = T)
      
      for (m in 1:mod_par$n_species) { # recording data in species data.frame
        output_species[l, m] <- sum(agents_grid@data@values == m, na.rm = T)
      }
      
      output$shannon[l] <- diversity(output_species[l, ], index = "shannon")
      
      clump_freq <- na.omit(freq(clumped)) # creating a frequency table from the patches raster
      
      n_spec_patch <- sapply(clump_freq[, 1], function(id) {
        length(unique(na.omit(agents_grid[clumped == id])))
      }) # creating a vector with number of species for each patch
      abundance_patch <- sapply(clump_freq[, 1], function(id) {
        length(na.omit(agents_grid[clumped == id]))
      }) # creating a vector with abundance value for each patch
      
      
      output_patches <- data.frame(
        sim_id = rep(counter, nrow(clump_freq)), # creating the patches data frame and adding the relevant data
        patch_id = clump_freq[, 1],
        area = clump_freq[, 2],
        n_species = n_spec_patch,
        individuals = abundance_patch
      )
      output_patches_all <- rbind(output_patches_all, output_patches) # binding data from output_patches to output_patches_all
    }
    output_all <- rbind(output_all, output) # binding data from output to output_all
    output_species_all <- rbind(output_species_all, output_species) # binding data from output_species ti output_species_all
  }
  
  output_all$sim_id <- seq(1, runs * nrow(var_par)) # assigning simulation ID to all rows in the general output file
  
  output_list <- list(general = output_all, species = output_species_all, patches = output_patches_all) # collecting all 3 output data.frame into a list
  
  end.time <- Sys.time() # for displaying run time of the function
  
  time.taken <- round(end.time - start.time, 2) # for displaying run time of the function
  print(time.taken)
  
  return(output_list)
}
