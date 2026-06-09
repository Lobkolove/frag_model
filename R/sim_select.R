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
  file_type = c("sampled", "state"),
  sampled = c("random", "checkerboard", "all"),
  log_file = "output/simulations_log.csv"
) {

  ftype <- match.arg(file_type)
  sampled <- match.arg(sampled)

  if (!file.exists(log_file)) stop("Log not found.")
  log <- fread(log_file)

  conds <- list(...)
  unused <- setdiff(names(conds), names(log))
  if (length(unused) > 0) stop("Unused conditions: ", paste(unused, collapse = ", "))

  matches <- log
  for (cond in names(conds)) {
    matches <- matches[matches[[cond]] == conds[[cond]], ]
  }
  if (nrow(matches) == 0) stop("No matches found.")
  
  if (ftype == "state") return(matches$state_file)
  if (ftype == "sampled") {
    pattern <- switch(sampled,
      "all" = "_all",
      "random" = "_rand",
      "checkerboard" = "_cb"
    )
    paths <- matches$sampled_file
    selected <- paths %>% 
      str_split("; ") %>% 
      unlist() %>% 
      str_subset(pattern = pattern)
    return(selected)
  }
}


sim_conditions <- function(log_file = "output/simulations_log.csv") {
  if (!file.exists(log_file)) stop("Log not found.")
  log <- fread(log_file)

  cat("Available conditions in log:\n")
  names <- setdiff(names(log), c("state_file", "sampled_files"))
  types <- sapply(log[, ..names], class)
  cat(paste0("- ", names, " (", types, ")\n"))
}
