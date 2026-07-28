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
    dispersal_dist = "disp"
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

  if (meta$dispersal_type == "random") {
    key_values[length(key_values)] <- "rand_disp"
  } 

  # Combine the key-value pairs into a single string.
  scenario_key <- paste0(key_values, collapse = "_")
  return(scenario_key)

}

rep_number <- function(
  meta,
  parameters = c("ac_amount", "habitat", "fragmentation", "edge_effect", "dispersal_type", "dispersal_dist"),
  log_file = "output/simulations_log.csv",
  log_dt = NULL
) {

  if (is.null(log_dt)) {
    if (!file.exists(log_file)) return(1L)
    log <- data.table::fread(log_file, colClasses = list(integer = c("sim_id", "replicate_num")))
  } else {
    log <- data.table::copy(log_dt)
  }
  
  if (nrow(log) == 0L) return(1L)
  if (!all(parameters %in% names(meta))) {
    stop("One or more specified parameters are not present in the metadata.")
  }

  if (is.null(log$dispersal_type)) {
    log$dispersal_type <- "short_long"
  }
  
  missing_cols <- setdiff(parameters, names(log))
  if (length(missing_cols) > 0) {
    stop("The following parameter columns are missing from the log file: ", paste(missing_cols, collapse = ", "))
  }

  # Build logical mask comparing each requested parameter
  mask <- rep(TRUE, nrow(log))
  for (par in parameters) {
    mask <- mask & (log[[par]] == meta[[par]])
  }
  current_id <- as.integer(meta$sim_id)
  mask <- mask & (log[["sim_id"]] < current_id)
  matches <- log[mask, ]

  if (nrow(matches) == 0) return(1L)
  return(max(matches$replicate_num, na.rm = TRUE) + 1L)

}

unique_sim_id <- function(
  log_file, 
  id_col = "sim_id",
  increment = 1L,
  as_string = FALSE, 
  format = "%04d"
) {

  if (!is.integer(increment)) stop("Argument 'increment' must be an integer")
  if (!file.exists(log_file)) {
    sim_id <- increment
  } else {
    log <- fread(log_file)
    last_sim_id <- max(as.numeric(log[[id_col]]))
    sim_id <- last_sim_id + increment
  }

  if (as_string) {
    return(sprintf(fmt = format, sim_id))
  } else {
    return(sim_id)
  }
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
  overwrite = FALSE,
  log_file = "output/simulations_log.csv"
) {
  current_id <- as.integer(meta$sim_id)

  if (file.exists(log_file)) {
    log <- fread(log_file, colClasses = list(character = c("sim_id", "replicate_num")))    
    if (current_id %in% as.integer(log$sim_id)) {
      if (!overwrite) {
        warning("Simulation ", current_id, " is already in the log and overwriting is disabled.\nSkipping log entry.\n", call. = FALSE)
        return()
      } else {
        warning("Simulation ", current_id, " is already in the log. Overwriting entry.\n", call. = FALSE)
        log <- log[as.integer(log$sim_id) != current_id]  # Remove existing entry for this sim_id
      }
    }
  } else {
    log <- data.table()  # Create an empty data.table if log file doesn't exist
  }

  # Check if meta contains run_date, if not assess from state file creation date
  if (!"run_date" %in% names(meta)) {
    run_date <- as.Date(file.info(state_file)$ctime)
  } else {
    run_date <- meta$run_date
  }

  entry <- data.table(
    sim_id = sprintf("%04d", current_id),
    job_id = job_id, 
    scenario_key = scenario_key, 
    replicate_num = sprintf("%03d", as.integer(replicate_num)),
    run_date = run_date, 
    project_version = project_version,
    master_seed = meta$master_seed, 
    ac_amount = meta$ac_amount, 
    fragmentation = meta$fragmentation, 
    habitat = meta$habitat, 
    niche_breadth = meta$niche_breadth, #Simon was here!
    dispersal_type = meta$dispersal_type, 
    dispersal_kernel = ifelse(meta$dispersal_type == "random", NA_character_, meta$dispersal_kernel),
    dispersal_ratio = ifelse(meta$dispersal_type == "random", NA_real_, meta$dispersal_ratio),
    dispersal_dist = ifelse(meta$dispersal_type == "random", NA_integer_, meta$dispersal_dist), 
    edge_effect = meta$edge_effect,
    status = status,
    state_file = as.character(fs::path_rel(state_file, start = here::here())),
    sampled_files = paste(fs::path_rel(sampled_files, start = here::here()), collapse = "; ")
  )

  # Insert entry into log file, keeping sorted by sim_id
  log <- rbind(log, entry)
  log <- log[order(as.integer(log$sim_id)), ]

  # Write updated log back to file
  fwrite(log, log_file)
  # fwrite(entry, log_file, append = file.exists(log_file))
  cat("\nA new entry for simulation", current_id, "was written to", log_file, "\n\n")
}
