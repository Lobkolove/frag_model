# First outline of a helper to query for specific simulations / scenarios
# This function is not working yet
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
  matches <- log %>% 
    dplyr::filter(conds)

  return(matches)

}
