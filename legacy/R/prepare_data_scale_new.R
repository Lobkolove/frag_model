require(dplyr)
require(tidyr)
prepare_data_scale_new <- function(path) {
  
  # get a list of all files matching a pattern
  
  data_all <- list.files(
    path = path,
    pattern = "output_general", full.names = TRUE
  )
  
  sim_info <- first(list.files(
    path = path,
    pattern = "static_parameters", full.names = TRUE
  ))
  sim_info_first <- read.csv(sim_info)

  #read the data
  
  read_data <- lapply(data_all, read.csv)
  data_1 <- bind_rows(read_data)
  
  # rename columns
  
  data <- dplyr::rename(data_1, abundance_landscape = individuals,
                        abundance_sample = sample_individuals_mean,
                        richness_landscape = present_species,
                        richness_sample = sample_species_mean,
                        diversity_landscape = shannon,
                        diversity_sample = sample_shannon_mean)

  # extract steps right after fragmentation and at the end of the simulation
  
  steps <- c(sim_info_first$V1[sim_info_first$X == "steps_pre_frag"] + 1, sim_info_first$V1[sim_info_first$X == "steps_post_frag"] + sim_info_first$V1[sim_info_first$X == "steps_pre_frag"])
  
  # extract data at the point right after fragmentation and at the end of simulation
  
  geo_data <- data[data$step == steps[1],]
  geo_demo_data <- data[data$step == steps[2],]
  
  # reshape data for plotting
  
  geo_rich <- tidyr::pivot_longer(geo_data, cols = c("richness_landscape", "richness_sample"), names_to = "richness_scale", values_to = "richness_val")
  gedo_rich <- tidyr::pivot_longer(geo_demo_data, cols = c("richness_landscape", "richness_sample"), names_to = "richness_scale", values_to = "richness_val")

  geo_rich$richness_scale <- paste0("geo_", geo_rich$richness_scale)
  gedo_rich$richness_scale <- paste0("gedo_", gedo_rich$richness_scale)
  combined_rich <- bind_rows(geo_rich, gedo_rich)
  
  landscape <- combined_rich %>%
    filter(grepl('landscape', richness_scale))
  
  sample <- combined_rich %>%
    filter(grepl("sample", richness_scale))
  
  list_of_outputs <- list(geo_rich = geo_rich, 
                          gedo_rich = gedo_rich, 
                          landscape = landscape,
                          sample = sample,
                          all = data_1)
  
  return(list_of_outputs)
}
