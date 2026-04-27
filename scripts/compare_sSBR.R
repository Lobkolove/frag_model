library(tidyverse)
library(scam)
source("R/sSBR.R")
source("R/ref_curve.R")

# Read in all sampled data files into a list
sampled_files <- list.files("output/sampled_data/", 
                            pattern = "samp_all.csv", 
                            full.names = TRUE)

data_list <- map(sampled_files, read_csv)
acc1 <- spat_acc(data_list[[1]], compute_richness = TRUE)

# Check ranges of spatial extent
ranges <- map_dfr(data_list, ~ {
  out <- spat_acc(.x, compute_richness = FALSE)
  tibble(
    fragmentation = unique(.x$fragmentation),
    max_spat_ext = max(out$spat_ext),
    n_samples = nrow(out)
  )
})
ranges

# Reformat each dataset to wide format
data_wide <- data_list %>%
  map(~ .x %>%
    dplyr::filter(step_label == "post_fragmentation") %>%
    tidyr::pivot_wider(names_from  = species_id,
                       values_from = n,
                       values_fill = 0,
                       names_prefix = "sp_",
                       names_sort = TRUE)
    )

# Build spatvec from 0 to maximum common area across all datasets (if not provided to ref_curve)


# Build reference curve using the median (0.5 quantile) of the area-based effort across all samples
eff_ref <- ref_curve(data_wide, 
                     quantile = 0.5, 
                     spatpar = "area", 
                     spatvec = NULL)

# # Merge together data which share the same fragmentation (data_list[[i]]$fragmentation)
# frag_levels <- unique(map_dbl(data_list, ~ .x$fragmentation[1]))
# data_merged <- vector("list", length(frag_levels))
# for (i in seq_along(frag_levels)) {
#   frag_level <- frag_levels[i]
#   # Bind first, pivot second to avoid NAs in species columns
#   data_merged[[i]] <- bind_rows(map(data_list, ~ filter(.x, fragmentation == frag_level))) %>%
#     tidyr::pivot_wider(names_from  = species_id,
#                        values_from = n,
#                        values_fill = 0,
#                        names_prefix = "sp_",
#                        names_sort = TRUE) %>%
#     select(-sp_NA)
# }

# Run sSBR for each dataset
out_list <- vector("list", length(data_wide))
for (i in seq_along(data_wide)) {
  start <- Sys.time()
  cat("Processing dataset", i, "( fragmentation =", data_wide[[i]]$fragmentation[1], "| ac_amount = " , data_wide[[i]]$ac_amount[1], ")\n")
  result <- sSBR(
    data_wide[[i]],
    method = "area",
    cutoff = FALSE,
    effort_ref = eff_ref
  )
  out_list[[i]] <- result$smooth
  end <- Sys.time()
  cat("Time taken: ", round(end - start, 2), "\n")
}

# Add fragmentation and autocorrelation metadata to each result
out_list <- map2(out_list, data_wide, ~ .x %>%
  mutate(
    fragmentation = .y$fragmentation[1],
    ac_amount = .y$ac_amount[1]
  )
)

# Combine results per autocorrelation level (3 results per level)
ac_amounts <- unique(map_dbl(data_wide, ~ .x$ac_amount[1]))
out_combined <- vector("list", length(ac_amounts))

for (i in seq_along(out_combined)) {
  ac_level <- ac_amounts[i]
  out_combined[[i]] <- bind_rows(map2(out_list, data_wide, ~ if (.y$ac_amount[1] == ac_level) .x else NULL)) %>%
    filter(!is.na(spat_ext))
}

# Plot results per autocorrelation level, color by fragmentation level
for (i in seq_along(out_combined)) {
  ac_level <- ac_amounts[i]
  p <- ggplot(out_combined[[i]], aes(x = spat_ext, y = S, color = factor(fragmentation))) +
    geom_line() +
    geom_ribbon(aes(ymin = S_low, ymax = S_high), alpha = 0.2) +
    labs(title = paste("Autocorrelation Level:", ac_level),
         x = "Spatial Extent (area)",
         y = "Predicted Richness (S)",
         color = "Fragmentation Level") +
    theme_minimal()
  
  print(p)
}

