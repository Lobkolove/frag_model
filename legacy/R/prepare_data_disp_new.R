require(broom)
require(tidyverse)

prepare_data_disp_new <- function(path) {
  
  # list and read data
  
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
  
  data_sample <- list.files(
    path = path,
    pattern = "output_sample", full.names = TRUE
  )
  
  var_par <- read.csv(var_par_read)
  
  sim_info_all <- read.csv(sim_info)
  sim_info_all <- as.data.frame(t(sim_info_all))
  
  read_data <- lapply(data_all, read.csv)
  data_1 <- bind_rows(read_data)
  
  read_data2 <- lapply(data_sample, read.csv)
  data_sam <- bind_rows(read_data2)
  
  # reshape
  
  colnames(sim_info_all) <- as.character(sim_info_all[1, ])
  sim_info_all <- sim_info_all[-1, ]
  colnames(var_par)[1] <- "sim_id"
  
  # extract data at the point right after fragmentation and at the end of simulation
  
  steps <- c(as.integer(sim_info_all$steps_pre_frag) + 1, as.integer(sim_info_all$steps_post_frag) + as.integer(sim_info_all$steps_pre_frag))
  
  data_1 <- data_1 %>% filter(step == steps[2] | step == steps[1])
  
  data_sel <- data_1 %>% 
    select(sim_id, step, fragmentation, habitat, present_species, sample_species_mean)
  
  var_par_dist <- var_par %>% select(sim_id, disp_dist)

  data <- left_join(data_sel, var_par_dist, by = "sim_id")
  
  # extract slopes
  
  data_step1 <- data %>% filter(step == steps[1])
  data_step2 <- data %>% filter(step == steps[2])
  
  data_lm_geo <- data_step1 %>%
    group_by(disp_dist, step) %>%
    nest() %>%
    mutate(lm = map(data, ~ lm(present_species ~ fragmentation, data = .x) %>%
                      tidy())) %>%
    unnest(lm) %>%
    filter(term == "fragmentation") %>%
    select(-term)
  
  data_lm_gedo <- data_step2 %>%
    group_by(disp_dist, step) %>%
    nest() %>%
    mutate(lm = map(data, ~ lm(present_species ~ fragmentation, data = .x) %>%
                      tidy())) %>%
    unnest(lm) %>%
    filter(term == "fragmentation") %>%
    select(-term)
  
  data_lm <- bind_rows(data_lm_gedo, data_lm_geo)
  
  output_list <- list(data = data,
                      data_lm = data_lm)
  
  return(output_list)
  
}