# create a data frame for species abundance and diversity results in the whole landscape
intialize_output_all <- function() {
  data.frame(
    sim_id = numeric(),
    step = numeric(),
    individuals = numeric(),
    present_species = numeric(),
    shannon = numeric(),
    patch_individuals_mean = numeric(),
    patch_species_mean = numeric(),
    patch_shannon_mean = numeric(),
    sample_individuals_mean = numeric(),
    sample_species_mean = numeric(),
    sample_shannon_mean = numeric(),
    habitat = numeric(),
    fragmentation = numeric(),
    n_patches = numeric(),
    n_samples = numeric()
  )
}

# create a data frame for species abundance and diversity results in the whole landscape
initialize_output_sample <- function(species_seq) {
  output_sample <- data.frame(
    sim_id = NA,
    patch_id = NA,
    sample_id = NA,
    step = NA,
    loc_x = NA,
    loc_y = NA,
    patch_size = NA,
    fragmentation = NA,
    habitat = NA
  )
  species_names <- paste0("sp_", species_seq)
  for (p in 1:length(species_names)) {
    output_sample[, as.character(species_names[p])] <- numeric()
  }
  return(output_sample)
}