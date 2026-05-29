# Helpers for filenaming and log management


scenario_key <- function(
  meta,
  parameters = c("ac_amount", "habitat", "fragmentation", "edge_effect", "dispersal_dist")
) {

  par_keys <- list(
    ac_amount = "ac",
    habitat = "hab",
    fragmentation = "frag",
    edge_effect = "edge",
    dispersal_dist = "disp",
    niche_breadth = "nb"
  )

  # Check if the specified parameters are valid and present in the metadata
  if (!all(parameters %in% names(meta))) {
    stop("One or more specified parameters are not present in the metadata.")
  }
  if (!all(parameters %in% names(par_keys))) {
    stop("One or more specified parameters do not have corresponding keys defined.\n\nParameters currently supported for the scenario key:\n", paste(names(par_keys), collapse = ", "))
  }

  # For the chosen parameters, extract their values from the metadata and create a key-value pair string.
  key_values <- sapply(parameters, function(par) {
    par_key <- par_keys[[par]]
    par_value <- meta[[par]]
    paste0(par_key, par_value)
  })

  # Combine the key-value pairs into a single string.
  scenario_key <- paste0(key_values, collapse = "_")
  return(scenario_key)

}


sim_filename <- function(sim_id, scenario_key, replicate_num) {
  if (is.numeric(sim_id)) {
    sim_id <- sprintf("%04d", sim_id)  # Pad with zeros to ensure 4 digits
  }
  if (is.numeric(replicate_num)) {
    replicate_num <- sprintf("%03d", replicate_num)  # Pad with zeros to ensure 3 digits
  }
  filename <- paste0(sim_id, "_", scenario_key, "_r", replicate_num)
  return(filename)
}


log_entry <- function(
  sim_id,
  job_id,
  scenario_key,
  replicate_num,
  project_version,
  master_seed,
  var_par,
  status,
  state_file,
  sampled_files,
  log_file = "output/simulations_log.csv"
) {

  entry <- data.table(
    sim_id = ifelse(is.numeric(sim_id), sprintf("%04d", sim_id), sim_id),
    job_id = job_id, 
    scenario_key = scenario_key, 
    replicate_num = ifelse(is.numeric(replicate_num), sprintf("%03d", replicate_num), replicate_num),
    run_date = Sys.Date(), 
    project_version = project_version,
    master_seed = master_seed, 
    ac_amount = var_par$ac, 
    fragmentation = var_par$frag, 
    habitat = var_par$hab, 
    niche_breadth = var_par$nb,
    dispersal = var_par$disp, 
    dispersal_dist = var_par$disp_dist, 
    edge_effect = var_par$edge,
    status = status,
    state_file = path_rel(state_file, start = here()),
    sampled_files = paste(path_rel(sampled_files, start = here()), collapse = "; ")
  )

  fwrite(entry, log_file, append = file.exists(log_file), logical01 = TRUE)
  cat("\nA log entry for simulation", sim_id, "was written to", log_file, "\n\n")
}
