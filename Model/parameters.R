################################################################################
#
# Setting model parameters
#
################################################################################

# choose a random seed constant for the simulation or set the seed manually
# seed <- round(runif(1,1,9999))
seed <- 2027

# Switches

switch <- data.frame(
  animation_export = 0, # determine whether plots from each time-step are exported for later animation. 0 = no, 1 = yes (keep OFF when running on cluster!)
  random = 0, # determine whether set.seed is used for the spin up phase. 0 = random, 1 = deterministic
  random_post_frag = 0, # determine whether set.seed is used for fragmenting landscape. 0 = random, 1 = non-random
  random_community = 0, # determine whether to introduce randomness in species distribution. 0 = random, 1 = non-random
  sample_all = 0, # number of samples according to sample parameter or sample all possible cells. 0 = parameter, 1 = all
  immigration = 1, # determine if individuals from outside the space can immigrate back in. 0 = no, 1 = yes

  ############
  # this switch is very problematic as it can override "var_par" in certain cases e.g with edge effects in "death" function.
  # This need to be examined and maybe removed. For now, leave at 0!!
  species_specific_par = 0, # determine if species have non/some/"all" unique parameters 0 = deterministic, 1 = some, 2 = "all"
  ############

  kernel_type = 1, # choose type of dispersal kernel.  0 = log-normal distribution, 1 = Exponential Distribution
  edge_effect = 0, # 0 = no edge effects 1 = with edge effects
  print_agents = 1, # for debugging. if switch = 1 a message with amount of agents in each step is printed
  export_raster = 0 # 0 = don't export gri and grd files of samples locations, 1 = export
)

if (switch$random == 1) {
  set.seed(seed)
}

# Static Parameters

mod_par <- data.frame(
  grid_size = 50, # side length of a square grid
  habitat_percent = 0.15, # 0-1 proportion of habitat vs matrix
  spatial_ac = 0.9, # autocorrelation of habitat 0 (rough) - 1 (smooth)
  frag_factor = 0.7, # level of fragmentation 0 (unified) - 1 (fragmented)
  n_pop = 5000, # setting amount of individuals
  n_species = 1000, # setting number of species
  niche_breadth = 0.1, # is used to determine the SD (nb) in exp((-(e-u)^2)/(2*nb^2))
  steps = 40, # Determine how many time steps in the dynamic model
  steps_pre_frag = 40, # used in GeDo_run.R instead of 'steps'
  steps_post_frag = 0, # used in GeDo_run.R instead of 'steps'
  birth_rate = 0.85, # chances of an individual giving birth
  death_rate = 0.25, # chances of an individual dying
  dispersal = 1, # Determining Long and short dispersal. value is proportion of short dispersal (0-1)
  mean_disp = 2, # Parameter for the mean dispersal distance
  sd_disp = 2, # parameter for the standard deviation for dispersal kernel
  k_inter = 50, # cell carrying capacity for all species
  k_intra = 50, # cell carrying capacity of same species individuals
  n_immigrants = 50, # number of immigrants per time-step. At the moment static, consider changing it
  n_samples = 30, # How many samples to collect
  pos_neg_edge = 1 # * with death rate, so EE < 1 is positive (reduce DR) and EE > 1 is negative
)

n_values <- seq(from = 0, to = 1, length.out = mod_par$n_species)

# creating species specific parameters instead of the static ones. This can be
# overridden with a switch

if (switch$species_specific_par == 2) {
  species_par <- data.frame(
    species_id = 1:mod_par$n_species,
    n_value = sample(n_values),
    birth_rate = rnorm(mod_par$n_species, mean = mod_par$birth_rate, sd = 0.1),
    death_rate = rnorm(mod_par$n_species, mean = mod_par$death_rate, sd = 0.1),
    dispersal_rate = rnorm(mod_par$n_species, mean = mod_par$dispersal, sd = 0.1),
    niche_breadth = rnorm(mod_par$n_species, mean = mod_par$niche_breadth, sd = 0.1),
    edge_effect = rnorm(mod_par$n_species, mean = mod_par$pos_neg_edge, sd = 0.1)
  )
} else if (switch$species_specific_par == 1) {
  species_par <- data.frame(
    species_id = 1:mod_par$n_species,
    n_value = sample(n_values),
    birth_rate = rep(mod_par$birth_rate, mod_par$n_species),
    death_rate = rep(mod_par$death_rate, mod_par$n_species),
    dispersal_rate = rep(mod_par$dispersal, mod_par$n_species),
    niche_breadth = rep(mod_par$niche_breadth, mod_par$n_species),
    edge_effect = rep(mod_par$pos_neg_edge, mod_par$n_species)
  )
} else if (switch$species_specific_par == 0) {
  species_par <- data.frame(
    species_id = 1:mod_par$n_species,
    n_value = n_values,
    birth_rate = rep(mod_par$birth_rate, mod_par$n_species),
    death_rate = rep(mod_par$death_rate, mod_par$n_species),
    dispersal_rate = rep(mod_par$dispersal, mod_par$n_species),
    niche_breadth = rep(mod_par$niche_breadth, mod_par$n_species),
    edge_effect = rep(mod_par$pos_neg_edge, mod_par$n_species)
  )
} else {
  print("check switches please")
}

# switches for var par to determine which parameters are constant and which aren't.
# if switch is on (1) the sequence of values will be used. off (0) will mean a constant value

vp_switch <- data.frame(
  frag = 0,
  ac = 0,
  hab = 0,
  nb = 0,
  disp = 0,
  disp_dist = 0,
  edge = 0
)
# Creating a table of varying parameter values using expand.grid() for the multiple runs function.
# Only fragmentation, auto-correlation, and habitat percent are adjustable and will override the settings of 'par'.
# To vary other parameters the 'multi_runs' function will have to be adjusted.

if (vp_switch$frag == 1) {
  frag_factor_vector <- seq(0.1, 0.9, 0.1)
} else {
  frag_factor_vector <- mod_par$frag_factor
}
if (vp_switch$ac == 1) {
  spatial_ac_vector <- seq(0.01, 0.91, 0.10)
} else {
  spatial_ac_vector <- mod_par$spatial_ac
}
if (vp_switch$hab == 1) {
  habitat_percent_vector <- seq(0.05, 0.95, 0.1)
} else {
  habitat_percent_vector <- mod_par$habitat_percent
}
if (vp_switch$nb == 1) {
  NB_vector <- seq(0.01, 0.3, 0.03)
} else {
  NB_vector <- mod_par$niche_breadth
}
if (vp_switch$disp == 1) {
  disp_vector <- disp_vector <- seq(0, 0.9, 0.1)
} else {
  disp_vector <- mod_par$dispersal
}
if (vp_switch$disp_dist == 1) {
  disp_dist_vector <- seq(1, 5.5, 0.5)
} else {
  disp_dist_vector <- mod_par$mean_disp
}
if (vp_switch$edge == 1) {
  edge_effect_vector <- seq(0.6, 1.5, 0.1)
} else {
  edge_effect_vector <- mod_par$pos_neg_edge
}

# Creates a table of all possible parameter combinations
var_par <- expand.grid(
  frag = frag_factor_vector,
  ac = spatial_ac_vector,
  hab = habitat_percent_vector,
  nb = NB_vector,
  disp = disp_vector,
  disp_dist = disp_dist_vector,
  edge = edge_effect_vector
)
