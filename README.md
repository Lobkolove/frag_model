# Disentangling geometric and demographic effects of habitat fragmentation using distance decay of similarity

*This repository contains the code used for my thesis and final presentation on the use of distance decay of similarity for disentangling geometric and demographic effects of habitat fragmentation from snapshot data. Below, you'll find details on how to run the metacommunity simulations and how to reproduce the diversity and distance decay analyses in
`scripts/main_analyses.R`.*


## Abstract

Habitat loss is among the leading causes of declines in terrestrial biodiversity, and what remains is typically left in smaller, more isolated patches. Whether that fragmentation is itself harmful remains disputed, in part because its effects depend on the scale at which diversity is measured and on the time elapsed since habitat was removed. Much of that disagreement can be resolved by separating geometric effects, which follow from spatial configuration alone, from demographic effects, which accumulate as ecological processes act on the fragmented landscape. Drawing the distinction, however, requires observing a landscape immediately after modification and again once those processes have acted, which empirical studies can rarely do. We show that the distance decay of similarity can make the same separation from a single point in time, because comparing similarity at matched distances leaves the curve insensitive to geometric effects. Using a spatially explicit, individual based model, we simulated communities on landscapes which differed in configuration at a constant amount of habitat, and analysed distance decay alongside alpha and gamma diversity. Curves from landscapes differing in configuration coincided immediately after habitat removal, although landscape species richness differed substantially between them at that moment. By the end of the simulations they had diverged, with more fragmented landscapes showing higher similarity between nearby sites and a steeper decline with distance. The pattern was recovered from samples of 30 sites, though with considerable variation between replicates, and its magnitude depended on how much habitat remained, on how strongly the environment was autocorrelated and on how far organisms dispersed. Distance decay can therefore isolate demographic effects in landscapes which can only be observed once. Used alongside local and landscape diversity, it opens a route to testing fragmentation effects with the data field studies already collect.

## Structure of the repository

```
Model/                  Simulation machinery
  parameters.R            Parameter grid and switches
  sim_array.sh            SLURM batch script (cluster entry point)
  sim_array.R             One simulation per array task (for cluster)
  sim_loop.R              Local loop over the parameter grid
  last_sim_id.R           Reads the last used sim_id from the log
  src/                    Model processes (birth, death, dispersal, immigration, ...)
R/                      Functions
  sample_cells.R          Sampling of habitat cells 
  toroidal_clump.R        Patch delineation on a toroidal grid
  sim_select.R            Query the simulation log for output files
  compute_diversity.R     Alpha/gamma richness, Hill numbers, evenness
  dist_decay.R            Distance decay of similarity
  toroidal_dist.R         Toroidal pairwise distances
scripts/
  main_analyses.R         The diversity and distance-decay analyses
output/
  simulations_log.csv     One row per simulation (parameters, seeds, file paths)
  model_states/           Full state snapshots (.rds), one per simulation
  sampled_data/           Sampled community data (.csv), the analysis input
pics/                   Exported figures
```

`output/model_states/`, `output/sampled_data/` and `pics/` are in `.gitignore` —
the data and figures are not versioned and have to be transferred separately.

All paths are resolved with `here::here()`, anchored on the `.here` file in the
root directory. Keep that file.

### Required R packages

- **Model:** here, fs, data.table, raster, dplyr, tidyr, igraph, collapse,
  checkmate, withr, scales, tibble, rlang
  (`igraph` is not called directly but `raster::clump()` needs it)
- **Analysis:** here, data.table, dplyr, ggplot2, ggh4x, patchwork, stringr,
  tidyr, purrr, vegan, tibble, rlang, mgcv

---

# Part 1 — Running the simulations

## 1.1 What one simulation does

`single_run()` (in `Model/src/single_run.R`) runs a single simulation:

1. **Initialize** a `grid_size × grid_size` toroidal landscape with 100 % habitat,
   an autocorrelated environmental gradient (`ac`), and `n_pop` individuals of
   `n_species` distributed according to their niche optima.
2. **Pre-fragmentation phase** — `steps_pre_frag` steps. Each step is birth
   (offspring disperse and then compete for recruitment against the local and
   conspecific carrying capacities), death, and immigration, so the community
   reaches a quasi-equilibrium.
3. **Fragmentation event** — the landscape is reduced to `hab` habitat with
   fragmentation level `frag`; all individuals outside the remaining habitat are
   removed and the destroyed cells become matrix.
4. **Post-fragmentation phase** — `steps_post_frag` further steps.

States can be recorded at four points: `start`, `pre_fragmentation`,
`post_fragmentation` and `final`. The contrast between `post_fragmentation`
(geometric effect only) and `final` (geometric + demographic effects) is what the
analysis is built on.

`single_run()` returns a named list of state snapshots; each snapshot holds the
landscape `grid`, the `agents` table, per-cell species abundances `ss_abund`, and
a `meta` list with all parameters and the master seed.

## 1.2 Setting the parameters

Everything is configured in `Model/parameters.R`:

- **`switch`** — model behaviour. The relevant ones here:
  - `dispersal_type`: `1` = short/long distance dispersal, `0` = dispersal to a
    random habitat cell (used for the "random dispersal" series).
  - `kernel_type`: `1` = exponential kernel, `0` = log-normal.
  - `edge_effect`: `0` = off, `1` = on (kept on `0` throughout the main analyses).
  - `species_specific_par`: `0` = same parameters for all species, `1` = species specific parameters (kept on `0` throughout the main analyses).
- **`mod_par`** — static parameters (grid size, number of species and individuals,
  birth/death rates, carrying capacities, number of steps, and the default values
  for habitat amount, autocorrelation, niche breadth and dispersal distance).
- **`vp_switch`** — which parameters vary across simulations. Per parameter:
  `0` = fixed at the `mod_par` value, `1` = three levels (low/medium/high),
  `2` = a full sequence.
- **`var_par`** — the resulting parameter grid, built with `expand.grid()`.
  One row = one scenario. **`nrow(var_par)` is what the array size has to match.**

The default `master_seed` is also set here. It is used as such by `sim_loop.R`
(incremented per simulation) but overridden in `sim_array.R`, where master seeds are generated based on array job and task IDs.

## 1.3 Running on the cluster (`sim_array.sh` → `sim_array.R`)

This is the entry point used for all simulations in the analysis.

```bash
cd Model
sbatch sim_array.sh
```

Submit from inside `Model/` — `sim_array.sh` calls `last_sim_id.R` and
`sim_array.R` by relative path.

Before submitting, check:

- `#SBATCH --array=1-N` in `sim_array.sh`. **N must be a multiple of
  `nrow(var_par)`.** Each pass through the grid is one replicate of every
  scenario, so `--array=1-18` with 9 scenarios gives 2 replicates each.
- `--time`, `--mem-per-cpu` and `--mail-user`.
- The settings block at the top of `sim_array.R`:
  - `steps_to_record` — which snapshots to keep.
  - `export_states` — write the full state list as `.rds` (`TRUE` by default).
  - `sampling_methods` — `all` (every habitat cell), `cb` (checkerboard) and/or
    `rand` (30 random habitat cells). One CSV per method.
  - `to_log` — must be `TRUE` for anything the analysis will use.
  - `spacer` — set to the size of a concurrently running array to keep `sim_id`s
    from colliding.

How the bookkeeping works:

- `last_sim_id.R` reads the highest `sim_id` from the log before the array starts;
  each task takes `sim_id = last_sim_id + task_id + spacer`.
- The master seed is `job_id %% 1000 + task_id`, so seeds hardly ever repeat
  across arrays. **The seeds are the identity of a replicate — the analysis
  checks for duplicated seeds before using a series** (see 2.3).
- The scenario is written into the filename by `scenario_key()` and
  `sim_filename()`:
  `0001_ac0.7_hab0.15_frag0.2_edge1_nb0.1_disp2_r001_all.csv`
  = `sim_id` `_` scenario `_r` replicate `_` sampling method.
  Random-dispersal runs get `dispR` instead of `disp<N>`.
- `log_entry()` appends one row per simulation to `output/simulations_log.csv`,
  with all parameters, the master seed, and the paths of the state and sampled
  files. This log is the index the analysis queries — nothing else links the
  files to their parameters.

## 1.4 Running locally (`sim_loop.R`)

`Model/sim_loop.R` runs the same model over `var_par` in a plain loop, for
testing or for a single run (keep `var_par` at one row).

```r
source(here::here("Model/sim_loop.R"))
```

Its defaults are the testing defaults: `export_states = FALSE`,
`sampling_methods = list()` (nothing sampled), `to_log = FALSE`. Uncomment the
sampling methods and set `to_log = TRUE` if the runs are meant to be used later.
Note that with `to_log = FALSE` the `sim_id` is just the loop index, which will
collide with existing files.

## 1.5 The parameter settings behind the five analysed series

| Series | Varied | Settings | sim_ids |
|---|---|---|---|
| Core | fragmentation | `vp_switch$frag = 1`, all others `0`; `mod_par$spatial_ac = 0.7`, `habitat_percent = 0.15`, `mean_disp = 2` | 0001–0060 |
| Habitat amount | fragmentation × habitat | `vp_switch$frag = 1`, `hab = 1` (0.1/0.3/0.5) | 0088–0177 |
| Dispersal distance | fragmentation × dispersal distance | `vp_switch$frag = 1`, `disp_dist = 1` (1/4/8 cells) | 0061–0087, 0178–0249 |
| Autocorrelation | fragmentation × autocorrelation | `vp_switch$frag = 1`, `ac = 1` (0/0.5/1) | 0250–0321, 0340–0357 |
| Random dispersal | fragmentation | `vp_switch$frag = 1`; `switch$dispersal_type = 0` | 0385–0414 |

Each series has 3 fragmentation levels × 10–20 replicates. The analysis reads the
`_all` samples of all five series, plus the `_rand` samples of the core series.

---

# Part 2 — The analysis (`scripts/main_analyses.R`)

## 2.1 Input

The script reads `output/simulations_log.csv` and the sampled CSVs in
`output/sampled_data/` — 429 files in total, ~3 GB. It does **not** read the
`.rds` state files, and it never touches the checkerboard (`_cb`) samples.

Each sampled CSV is in long format: one row per sampled cell and species, with
the cell coordinates (`x_loc`, `y_loc`), abundance (`n`), patch id and size, the
time step (`step_label`), and all simulation parameters repeated as columns.

Only `post_fragmentation` and `final` are used.

## 2.2 Selecting simulations

Two helpers in `R/sim_select.R` do the selection, both querying the log. For example:

```r
core_ids <- sim_ids(habitat == 0.15, ac_amount == 0.7, dispersal_dist == 2)
paths    <- sim_select(sim_id %in% core_ids, sampled = "all")
```

- `sim_ids(...)` takes any filter expression over log columns and returns the
  matching `sim_id`s.
- `sim_select(...)` returns file paths. `vars` names the parameter that is allowed
  to vary (default `"fragmentation"`); for every other parameter it keeps only one
  value. With `mode = "biggest_sample"` (the default) it silently keeps the most
  frequent combination — use `mode = "user_defined"` to see the alternatives and
  pick. `sampled` selects `"all"`, `"random"` or `"checkerboard"`;
  `file_type = "state"` returns the `.rds` paths instead.
- `sim_vars()` lists the queryable columns.

## 2.3 Seed check

Before each series the script pulls the `master_seed`s of the selected
simulations and counts duplicates. Duplicated seeds mean the same landscape and
starting community were reused, so those runs are not independent replicates.
Four of the five series are clean; in the autocorrelation series 18 duplicates
are found and dropped with `ac_ids[!duplicated(ac_seeds)]`. **Keep this check
when adding new series.**

## 2.4 The two response measures

**Distance decay** — `grouped_ddecay()` in `R/dist_decay.R` splits a dataset by
`group_cols` (fragmentation × time step, plus the series parameter) and, per
group, calls `dist_decay()`: pivots the long data to a site × species matrix,
computes pairwise Bray–Curtis similarity (`1 - vegdist`) and pairwise spatial
distance, and fits `similarity ~ s(distance)` with `mgcv::gam()`. Predictions are
returned on a common `distvec` — always pass one explicitly (the script uses
`seq(0, 25, length.out = 200)`) so curves are comparable across groups.
`dist_type = "toroidal"` uses `toroidal_dist()` instead of Euclidean distance;
the analysis uses the Euclidean default.

**Diversity indices** — `compute_diversity()` in `R/compute_diversity.R` returns,
per time step and per scale, four indices in long format: `richness`,
`hill_shannon`, `hill_simpson` and `evenness`. *Local scale* is the mean over
sampled cells, *landscape scale* is the pooled community. Pass the series
parameter via `metadata_cols` (e.g. `metadata_cols = "habitat"`) so it is carried
into the result; it must have exactly one value per dataset.

Both are mapped over the list of datasets, `bind_rows()`ed, and summarised as mean
with 2.5 % / 97.5 % quantiles across replicates.

## 2.5 Running it

```r
source(here::here("scripts/main_analyses.R"))
```

The script is written to be stepped through section by section in RStudio rather
than sourced in one go — each section prints its plots and its seed check.
Sections are independent, so a single series can be re-run on its own after the
setup block (libraries, `source()` calls, palettes, log) at the top.

Expect it to take a while: 429 CSVs are read into memory and one GAM is fitted per
group.

Failures in `grouped_ddecay()` are caught per dataset in the dispersal,
autocorrelation and random series, and the failed index is reported. In the
dispersal series, dataset 54 fails and the whole replicate (datasets 46–54) is
dropped with `dd_dispersal[46:54] <- NULL`. **This index is positional — if the
set of input files changes, re-check which datasets fail instead of reusing the
hard-coded range.**

## 2.6 Output figures

Six figures are written to `pics/` (300 dpi PNG):

| File | Content |
|---|---|
| `core_combined.png` | Core series, all habitat cells: distance decay + all four indices |
| `core_random_combined.png` | Core series, 30 random cells: same layout |
| `diversity_core_random_combined.png` | Indices for full vs. random sample, stacked |
| `habitat_combined_aligned.png` | Habitat amount series: distance decay + richness |
| `dispersal_combined_aligned.png` | Dispersal distance series: distance decay + richness |
| `ac_combined_aligned.png` | Autocorrelation series: distance decay + richness |

Section 5 (random dispersal) only prints its plots; nothing is saved.

Plot conventions, set once at the top of the script: fragmentation is
low/medium/high in blue/yellow/purple (`pal_frag`), time step is
post-fragmentation vs. end of simulation in dark violet/pink (`pal_geodem`),
`theme_bw(base_size = 14)`, panels tagged A/B with `patchwork`.

---
