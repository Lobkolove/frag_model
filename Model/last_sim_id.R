#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop("Usage: last_sim_id.R LOG_FILE")
}

log_file <- args[[1]]

if (!file.exists(log_file) || file.info(log_file)$size == 0) {
  cat(0L)
  quit(status = 0L)
}

log <- read.csv(log_file, stringsAsFactors = FALSE)

if (!"sim_id" %in% names(log)) {
  stop("Column 'sim_id' not found in: ", log_file)
}

if (nrow(log) == 0L) {
  cat(0L)
} else {
  sim_ids <- as.integer(log$sim_id)
  sim_ids <- sim_ids[!is.na(sim_ids)]

  if (length(sim_ids) == 0L) {
    cat(0L)
  } else {
    cat(max(sim_ids))
  }
}