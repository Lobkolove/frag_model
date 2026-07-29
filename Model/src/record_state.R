record_step <- function(
  grid,
  agents,
  step,
  cells = NULL
) {
  structure(
    list(
      grid = grid,
      agents = agents,
      step = step,
      cells = cells
    ),
    class = c("model_state", "core_state", "list")
  )
}


build_meta <- function(
  meta = NULL,
  ...
) {

  if (!is.null(meta) && !is.list(meta)) stop("`meta` must be a list or NULL.")
  out <- if (is.null(meta)) list() else meta

  required <- c(
    "sim_id", "run_date", "master_seed", "grid_size",
    "ac_amount", "habitat", "fragmentation", "niche_breadth", "edge_effect", 
    "dispersal_type", "dispersal_kernel", "dispersal_ratio", "dispersal_dist"
  )

  provided <- list(...)

  if (length(provided) > 0L) {
    unused <- setdiff(names(provided)[names(provided) != "edge"], required)

    if (length(unused) > 0L) warning("One or more provided parameters were unexpected and excluded from the output meta object: ", unused)
    out[names(provided)] <- provided
  } 

  if (is.null(out$run_date)) out$run_date <- Sys.Date()
  if (is.null(out$edge_effect) && !is.null(out$edge)) out$edge_effect <- out$edge
  if (is.null(out$dispersal_type)) out$dispersal_type <- "short_long"
  if (is.null(out$dispersal_kernel)) out$dispersal_kernel <- "exponential"
  if (is.null(out$dispersal_ratio ) && !is.null(out$dispersal)) out$dispersal_ratio <- out$dispersal

  missing <- required[vapply(required, function(x) is.null(out[[x]]), logical(1))]
  if (length(missing) > 0L) stop("One or more required parameters are missing:\n", missing)

  out <- out[required]
  class(out) <- c("meta", "list")
  return(out)
}


# Add metadata and species abundance table to an existing core_state.
record_state <- function(
  model_state,
  meta = NULL,
  step_label
) {
  
  force(step_label)

  if (!is.list(model_state)) {
    stop("`model_state` must be a list, preferably a model state object.")
  }

  if (is.null(meta)) {
    if (is.null(model_state$meta)) stop(
    "No metadata found. If you're trying to record a core model state with no metadata, make sure to provide a metadata object using the `meta` argument."
    ) 
    meta <- model_state$meta
  }

  required <- c("grid", "agents", "step")

  missing <- setdiff(required, names(model_state))
  if (length(missing) > 0L) {
    stop("`model_state` missing element(s): ", paste(missing, collapse = ", "))
  }

  nulls <- required[vapply(required, function(x) is.null(x), logical(1))]
  if (length(nulls) > 0L) {
    stop("`model_state` has empty element(s): ", paste(nulls, collapse = ", "))
  }

  if (!is.null(model_state$ss_abund)) {
    ss_abund <- model_state$ss_abund
  } else {
    ss_abund <- dplyr::count(model_state$agents, x_loc, y_loc, species_id, name = "n")
  }

  full_state <- structure(
    list(
      meta = meta,
      step = model_state$step,
      step_label = step_label,
      grid = model_state$grid,
      agents = model_state$agents,
      cells = model_state$cells,
      ss_abund = ss_abund
    ),
    class = c("model_state", "full_state", "list")
  )

  return(full_state)
}
