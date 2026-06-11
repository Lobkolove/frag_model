#' Select simulation output files by matching conditions
#'
#' Queries a simulation log file to find outputs matching specified conditions,
#' then returns paths to either state files or sampled output files. For sampled
#' outputs, you can specify which sampling strategy (random, checkerboard, or all)
#' to retrieve.
#'
#' @param ... Condition pairs to match against log columns (e.g., `param1 = value1, param2 = value2`).
#'   All conditions must be valid column names in the log file. To get a list of valid conditions, see `sim_conditions()`.
#' @param file_type Type of output file to return. One of `"sampled"` (default) or `"state"`.
#' @param sampled For sampled files, which sampling strategy to return: `"random"`, `"checkerboard"`,
#'   or `"all"`. Ignored if `file_type = "state"`.
#' @param log_file Path to the simulation log CSV file. Defaults to `"output/simulations_log.csv"`.
#'
#' @return A character vector of file paths matching the specified conditions.
#'
#' @details
#' The function expects the log file to contain columns for each condition as well as
#' `state_file` and `sampled_file` columns. When `file_type = "sampled"`, the `sampled_file`
#' column should contain semicolon-separated paths (e.g., `"path_rand; path_cb; path_all"`).
#'
#' @examples
#' \dontrun{
#' # Get state files for a specific set of conditions
#' sim_select(fragmentation = 0.7, dispersal_dist = 10, file_type = "state")
#'
#' # Get randomly sampled output files
#' sim_select(fragmentation = 0.7, dispersal_dist = 10, file_type = "sampled", sampled = "random")
#' }
#'
#' @importFrom data.table fread
#' @importFrom stringr str_split str_subset
#' @importFrom magrittr %>%
#' @export
sim_select <- function(
  ...,
  vars = "fragmentation",
  file_type = c("sampled", "state"),
  sampled = c("random", "checkerboard", "all"),
  ignore = c(
    "sim_id",
    "job_id",
    "scenario_key",
    "run_date",
    "project_version",
    "master_seed",
    "replicate_num"
  ),
  mode = c("biggest_sample", "user_defined"),
  log_file = here("output/simulations_log.csv")
) {

  # Validate arguments
  ftype <- match.arg(file_type)
  sampled <- match.arg(sampled)
  mode <- match.arg(mode)

  # Check that log file exists and read it in
  if (!file.exists(log_file)) {
    stop("Log not found.")
  }
  log <- fread(log_file)

  # Check that specified variables and conditions are valid columns in the log
  unused <- setdiff(c(vars, ignore), names(log))
  if (length(unused) > 0) {
    stop(
      "Unused variables: ",
      paste(unused, collapse = ", "),
      "\n To see valid variables, use sim_vars()."
    )
  }

  # First, filter matches by the user-specified conditions.
  # If no conditions are specified, we will return all matches for now.
  matches <- log %>%
    filter(...) %>% 
    select(-all_of(ignore))
  if (nrow(matches) == 0) stop("No matches found.")
  
  # Now we want to filter matches so that the static variables are the same across all returned rows.
  # Whenever we encounter multiple values for a static variable, we will either:
  # mode == "biggest_sample": keep only the most common combination of all static variables;
  # mode == "user_defined": ask the user to specify the wanted combination. Group sizes should still 
  # be printed to help the user make an informed decision.
  static <- matches %>%
    select(-all_of(vars), -state_file, -sampled_files) %>%
    names()
  combs <- matches %>%
    group_by(across(all_of(static))) %>%
    summarise(n = n(), .groups = "drop")
  if (mode == "biggest_sample") {
    biggest <- combs %>%
      filter(n == max(n)) %>%
      select(-n)
    matches <- matches %>%
      inner_join(biggest, by = static)
  } else if (mode == "user_defined") {
    cat("Multiple combinations of static variables found:\n")
    print(combs)
    cat("\nPlease specify which combination to use by entering the row number:\n")
    user_input <- as.numeric(readline())
    if (is.na(user_input) || user_input < 1 || user_input > nrow(combs)) {
      stop("Invalid row number.")
    }
    matches <- matches %>%
      inner_join(select(combs[user_input,], -n), by = static)
  }


  # Finally, we can return the requested file paths. 
  # If file_type == "state", we will return the state_file column. 
  # If file_type == "sampled", we will return the sampled_file column, 
  # but only the paths that match the specified sampling strategy.
  if (ftype == "state") {
    return(matches$state_file)
  }
  if (ftype == "sampled") {
    pattern <- switch(
      sampled,
      "all" = "_all",
      "random" = "_rand",
      "checkerboard" = "_cb"
    )
    paths <- matches$sampled_files
    selected <- paths %>%
      str_split("; ") %>%
      unlist() %>%
      str_subset(pattern = pattern)
    return(selected)
  }
}


sim_vars <- function(log_file = "output/simulations_log.csv") {
  if (!file.exists(log_file)) stop("Log not found.")
  log <- fread(log_file)

  cat("Available variables in log:\n")
  names <- setdiff(names(log), c("state_file", "sampled_files"))
  types <- sapply(log[, ..names], class)
  cat(paste0("- ", names, " (", types, ")\n"))
}

sim_ids <- function(
  ..., 
  log_file = here("output/simulations_log.csv")
) {

  qs <- dplyr::enquos(...)
  used_vars <- unique(unlist(lapply(qs, function(q) all.vars(rlang::quo_get_expr(q)))))
  print(used_vars)
  
  # Check that log file exists and read it in
  if (!file.exists(log_file)) {
    stop("Log not found.")
  }
  log <- fread(log_file)

  matches <- log %>% 
    dplyr::filter(...)

  return(matches$sim_id)

}
