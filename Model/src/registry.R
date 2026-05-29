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

rep_number <- function(
  meta,
  parameters = c("ac_amount", "habitat", "fragmentation", "edge_effect", "dispersal_dist"),
  log_file = "output/simulations_log.csv"
) {

  if (!file.exists(log_file)) return(1L) 
  
  if (!all(parameters %in% names(meta))) {
    stop("One or more specified parameters are not present in the metadata.")
  }
    
  log <- data.table::fread(log_file, colClasses = list(integer = "replicate_num"))
  
  missing_cols <- setdiff(parameters, names(log))
  if (length(missing_cols) > 0) {
    stop("The following parameter columns are missing from the log file: ", paste(missing_cols, collapse = ", "))
  }

  # Build logical mask comparing each requested parameter
  mask <- rep(TRUE, nrow(log))
  for (par in parameters) {
    mask <- mask & (log[[par]] == meta[[par]])
  }
  matches <- log[mask, ]

  if (nrow(matches) == 0) return(1L)
  return(max(matches$replicate_num, na.rm = TRUE) + 1L)

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
  meta,
  job_id,
  scenario_key,
  replicate_num,
  project_version,
  status,
  state_file,
  sampled_files,
  log_file = "output/simulations_log.csv"
) {
  
  sim_id <- ifelse(is.numeric(meta$sim_id), sprintf("%04d", meta$sim_id), meta$sim_id)
  entry <- data.table(
    sim_id = sim_id,
    job_id = job_id, 
    scenario_key = scenario_key, 
    replicate_num = ifelse(is.numeric(replicate_num), sprintf("%03d", replicate_num), replicate_num),
    run_date = Sys.Date(), 
    project_version = project_version,
    master_seed = meta$master_seed, 
    ac_amount = meta$ac_amount, 
    fragmentation = meta$fragmentation, 
    habitat = meta$habitat, 
    niche_breadth = meta$niche_breadth, #I was here!
    dispersal = meta$dispersal, 
    dispersal_dist = meta$dispersal_dist, 
    edge_effect = meta$edge_effect,
    status = status,
    state_file = path_rel(state_file, start = here()),
    sampled_files = paste(path_rel(sampled_files, start = here()), collapse = "; ")
  )

  fwrite(entry, log_file, append = file.exists(log_file), logical01 = TRUE)
  cat("\nA log entry for simulation", sim_id, "was written to", log_file, "\n\n")
}
