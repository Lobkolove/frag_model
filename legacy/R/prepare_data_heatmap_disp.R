require(tidyverse)
require(broom)

prepare_data_heatmap_disp <- function(path, var1, var2) {
  # get a list of all files matching a pattern

  data_all <- list.files(
    path = path,
    pattern = "output_general", full.names = TRUE
  )

  sim_info <- list.files(
    path = path,
    pattern = "static_parameters", full.names = TRUE
  )

  var_par_read <- list.files(
    path = path,
    pattern = "varaying_parameters", full.names = TRUE
  )
  # read the data

  var_par <- read.csv(var_par_read)

  sim_info_all <- read.csv(sim_info)
  sim_info_all <- as.data.frame(t(sim_info_all))

  # reshape
  colnames(sim_info_all) <- as.character(sim_info_all[1, ])
  sim_info_all <- sim_info_all[-1, ]
  colnames(var_par)[1] <- "sim_id"

  read_data <- lapply(data_all, read.csv)
  data_1 <- bind_rows(read_data)

  # extract steps right after fragmentation and at the end of the simulation -needed?

  steps <- c(as.integer(sim_info_all$steps_pre_frag) + 1, as.integer(sim_info_all$steps_post_frag) + as.integer(sim_info_all$steps_pre_frag))



  ### prepare the data ###

  # extract the data at the 2 time points

  data_geo <- data_1 %>% filter(step == steps[1])
  data_geo_1 <- full_join(var_par, data_geo, by = "sim_id")
  data_geo_1 <- data_geo_1 %>%
    select(-n_patches, -n_samples, -individuals, -shannon, -patch_individuals_mean, -patch_species_mean, -patch_shannon_mean, -sample_individuals_mean, -sample_species_mean, -sample_shannon_mean) %>%
    rename(richness = present_species)
  data_gedo <- data_1 %>% filter(step == steps[2])
  data_gedo_1 <- full_join(var_par, data_gedo, by = "sim_id")
  data_gedo_1 <- data_gedo_1 %>%
    select(-n_patches, -n_samples, -individuals, -shannon, -patch_individuals_mean, -patch_species_mean, -patch_shannon_mean, -sample_individuals_mean, -sample_species_mean, -sample_shannon_mean) %>%
    rename(richness = present_species)
  
  # all is left to do - calculate slope for richness values in each simulation and add to table. Then group the table according to scenario

  data_lm_geo <- data_geo_1 %>%
    group_by((!!as.symbol(var1)), (!!as.symbol(var2))) %>%
    nest() %>%
    mutate(lm = map(data, ~ lm(richness ~ fragmentation, data = .x) %>%
      tidy())) %>%
    unnest(lm) %>%
    filter(term == "fragmentation") %>%
    select(-term)

  data_lm_gedo <- data_gedo_1 %>%
    group_by((!!as.symbol(var1)), (!!as.symbol(var2))) %>%
    nest() %>%
    mutate(lm = map(data, ~ lm(richness ~ fragmentation, data = .x) %>%
      tidy())) %>%
    unnest(lm) %>%
    filter(term == "fragmentation") %>%
    select(-term)

  results <- list(geo_lm = data_lm_geo, gedo_lm = data_lm_gedo, data = data_1, steps = steps)
  return(results)
}

