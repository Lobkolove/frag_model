library(data.table)
library(here)
library(stringr)
library(tools)

state_dir <- "output/model_states"
sampled_dir <- "output/sampled_data"
output_dir <- here("output")
log_path <- file.path(output_dir, "simulations_log.csv")

parse_filename <- function(fname) {
  base_name <- file_path_sans_ext(basename(fname))
  pattern <- "^(\\d{4})_(ac[0-9.]+)_(frag[0-9.]+)_(hab[0-9.]+)_r(\\d{3})(?:_samp_.*)?$"
  matches <- str_match(base_name, pattern)
  if (is.na(matches[1])) return(NULL)
  sim_id <- matches[2]
  ac <- as.numeric(str_remove(matches[3], "ac"))
  frag <- as.numeric(str_remove(matches[4], "frag"))
  hab <- as.numeric(str_remove(matches[5], "hab"))
  rep_num <- as.integer(str_remove(matches[6], "r"))
  data.table(sim_id = sim_id, ac = ac, frag = frag, hab = hab, rep_num = rep_num)
}

# list files (full paths)
state_files <- list.files(state_dir, "\\.rds$", full.names = TRUE)
sampled_files <- list.files(sampled_dir, "\\.csv$", full.names = TRUE)

safe_parse <- function(path) {
  parsed <- parse_filename(path)
  if (is.null(parsed)) return(NULL)
  parsed[, file_path := path]
  parsed
}

# build states table (keep full path in state_file)
state_list <- Filter(Negate(is.null), lapply(state_files, safe_parse))
state_dt <- if (length(state_list) > 0) {
  rbindlist(state_list, fill = TRUE)[
    , .(sim_id, ac, frag, hab, replicate_num = rep_num,
        state_file = file_path,
        run_date = as.Date(file.mtime(file_path)),
        nb = 0.1, disp = 1, disp_dist = 2, edge = 1,
        project_version = "frag_v1", job_id = "recovered",
        status = "complete", seed = NA_integer_)
  ]
} else data.table()

# build sampled table grouped by simulation (keep full paths)
sampled_list <- Filter(Negate(is.null), lapply(sampled_files, safe_parse))
sampled_dt <- if (length(sampled_list) > 0) rbindlist(sampled_list, fill = TRUE) else data.table()
sampled_grouped <- if (nrow(sampled_dt) > 0) {
  sampled_dt[, .(
    sampled_files = paste(unique(file_path), collapse = "; "),
    run_date = as.Date(max(file.mtime(file_path)))
  ), by = .(sim_id, ac, frag, hab, replicate_num = rep_num)]
} else data.table()

# outer join state + sampled on canonical keys
if (nrow(state_dt) == 0 && nrow(sampled_grouped) == 0) {
  message("No recoverable files found")
} else {
  merged <- merge(
    state_dt, sampled_grouped,
    by = c("sim_id", "ac", "frag", "hab", "replicate_num"),
    all = TRUE, suffixes = c(".state", ".samp")
  )

  # scenario_key, unified run_date, fill defaults where missing
  merged[, scenario_key := paste0("ac", sprintf("%.1f", ac), "_frag", frag, "_hab", sprintf("%.2f", hab))]
  merged[, run_date := as.Date(fifelse(!is.na(run_date.state), run_date.state, run_date.samp))]
  merged[is.na(nb), nb := 0.1]
  merged[is.na(disp), disp := 1]
  merged[is.na(disp_dist), disp_dist := 2]
  merged[is.na(edge), edge := 1]

  # infer task_id from sim_id suffix (last up to 3 digits) -- best-effort
  merged[, task_id := as.integer(substr(sim_id, pmax(1, nchar(sim_id) - 2), nchar(sim_id)))]

  # build final log with exact column order used in extended_pipeline.R
  log_entries <- merged[, .(
    sim_id = sim_id,
    task_id = task_id,
    job_id = fifelse(!is.na(job_id), job_id, "recovered"),
    scenario_key = scenario_key,
    replicate_num = as.integer(replicate_num),
    run_date = as.Date(run_date),
    project_version = "frag_v1",
    ac = ac,
    frag = frag,
    hab = hab,
    nb = nb,
    disp = disp,
    disp_dist = disp_dist,
    edge = edge,
    seed = as.integer(NA),      # unknown unless in .rds
    status = "complete",
    state_file = fifelse(!is.na(state_file), state_file, NA_character_),
    sampled_files = fifelse(!is.na(sampled_files), sampled_files, NA_character_)
  )]

  # ensure sim_id is zero-padded character like pipeline
  log_entries[, sim_id := as.character(sim_id)]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(log_entries, log_path)
  message("✅ Recovered ", nrow(log_entries), " simulations -> ", log_path)
}