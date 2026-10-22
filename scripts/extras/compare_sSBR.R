library(tidyverse)
library(scam)
source("R/sSBR.R")
source("R/ref_curve.R")

# Read in all full sampled datasets with an autocorrelation of 0.7
sampled_files <- list.files("output/sampled_data/", 
                            pattern = "ac0\\.7.*samp_all\\.csv$", 
                            full.names = TRUE)

# Create a list with all datasets
data_list <- map(sampled_files, read_csv)

# Select post-frag only and reformat to wide format
data_wide <- data_list %>%
  map(
    ~ .x %>%
      dplyr::filter(step_label == "post_fragmentation") %>%
      tidyr::pivot_wider(
        names_from = species_id,
        values_from = n,
        values_fill = 0,
        names_prefix = "sp_",
        names_sort = TRUE
      )
  )
  
# Check ranges of spatial extent
ranges <- map_dfr(data_wide, ~ {
  out <- spat_acc(.x, compute_richness = FALSE)
  tibble(
    sim_id = .x$sim_id[1],
    fragmentation = unique(.x$fragmentation),
    max_spat_ext = max(out$spat_ext),
    n_samples = nrow(out)
  )
})
ranges


# Test 1: mean spat_ext (per sampling effort) throughout all simulations ------------------------

data_acc <- map_dfr(data_wide, ~ spat_acc(.x, compute_richness = FALSE, average = TRUE))

eff_ref <- data_acc %>%
  group_by(samp_eff) %>%
  summarise(spat_ext = mean(spat_ext, na.rm = TRUE), .groups = "drop")

ggplot(data_acc, aes(x = samp_eff, y = spat_ext)) +
  geom_line(aes(group = sim_id), alpha = 0.3) +
  geom_line(data = eff_ref, color = "red", linewidth = 1) +
  labs(x = "Sampling Effort", y = "Spatial Extent (area)") +
  theme_minimal()

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
  out_list[[i]] <- result$smooth %>% 
    mutate(sim_id = data_wide[[i]]$sim_id[1],
           fragmentation = data_wide[[i]]$fragmentation[1])
  end <- Sys.time()
  cat("Time taken: ", round(end - start, 2), "\n")
}

out_df <- bind_rows(out_list)

sSBR_avg <- out_df %>%
  group_by(fragmentation, spat_ext, samp_eff) %>%
  summarise(
    S_mean = mean(S), 
    S_low = quantile(S, 0.025),
    S_high = quantile(S, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    fragmentation = factor(fragmentation, levels = c(0.2, 0.5, 0.8),
                            labels = c("low", "medium", "high"))
  )

pal <- c(colorspace::lighten("midnightblue", 0.2), 
         colorspace::lighten("violetred4", 0.2), 
         colorspace::lighten("darkgoldenrod", 0.2))

ggplot(sSBR_avg, aes(x = spat_ext, y = S_mean, color = fragmentation, fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), alpha = 0.4, color = NA) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = colorspace::lighten(pal, 0.2)) +
  labs(x = "Cumulative area of convex hull", y = "Cumulative species richness") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

# Visualize with spatial extent cutoff at 625
sSBR_cutoff <- out_df %>%
  filter(spat_ext <= 625) %>% 
  group_by(fragmentation, spat_ext, samp_eff) %>%
  summarise(
    S_mean = mean(S), 
    S_low = quantile(S, 0.025),
    S_high = quantile(S, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    fragmentation = factor(fragmentation, levels = c(0.2, 0.5, 0.8),
                            labels = c("low", "medium", "high"))
  )

ggplot(sSBR_cutoff, aes(x = spat_ext, y = S_mean, color = fragmentation, fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), alpha = 0.2, color = NA) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = "Cumulative area of convex hull", y = "Cumulative species richness") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")



# Test 2 - median spat_ext -----------------------------------------------

eff_ref_median <- data_acc %>%
  group_by(samp_eff) %>%
  summarise(spat_ext = median(spat_ext, na.rm = TRUE), .groups = "drop")

out_df_median <- map_dfr(
  data_wide,
  ~ sSBR(
    .x,
    method = "area",
    cutoff = TRUE,
    effort_ref = eff_ref_median
  )$smooth %>% 
    mutate(sim_id = .x$sim_id[1],
           fragmentation = .x$fragmentation[1])
)

out_list_median <- vector("list", length(data_wide))
for (i in seq_along(data_wide)) {
  start <- Sys.time()
  cat("Processing dataset", i, "( fragmentation =", data_wide[[i]]$fragmentation[1], "| ac_amount = " , data_wide[[i]]$ac_amount[1], ")\n")
  result <- sSBR(
    data_wide[[i]],
    method = "area",
    cutoff = TRUE,
    effort_ref = eff_ref_median
  )
  out_list_median[[i]] <- result$smooth %>% 
    mutate(sim_id = data_wide[[i]]$sim_id[1],
           fragmentation = data_wide[[i]]$fragmentation[1])
  end <- Sys.time()
  cat("Time taken: ", round(end - start, 2), "\n")
}

out_df_median <- bind_rows(out_list_median)

sSBR_avg_median <- out_df_median %>%
  group_by(fragmentation, spat_ext, samp_eff) %>%
  summarise(
    S_mean = mean(S), 
    S_low = quantile(S, 0.025),
    S_high = quantile(S, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    fragmentation = factor(fragmentation, levels = c(0.2, 0.5, 0.8),
                            labels = c("low", "medium", "high"))
  )

ggplot(sSBR_avg_median, aes(x = spat_ext, y = S_mean, color = fragmentation, fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), alpha = 0.4, color = NA) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = colorspace::lighten(pal, 0.2)) +
  labs(x = "Cumulative area of convex hull", y = "Cumulative species richness") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

# Again with cutoff at 625
sSBR_cutoff_median <- out_df_median %>%
  filter(spat_ext <= 625) %>% 
  group_by(fragmentation, spat_ext, samp_eff) %>%
  summarise(
    S_mean = mean(S), 
    S_low = quantile(S, 0.025),
    S_high = quantile(S, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    fragmentation = factor(fragmentation, levels = c(0.2, 0.5, 0.8),
                            labels = c("low", "medium", "high"))
  )

ggplot(sSBR_cutoff_median, aes(x = spat_ext, y = S_mean, color = fragmentation, fill = fragmentation)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = S_low, ymax = S_high), alpha = 0.2, color = NA) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(x = "Cumulative area of convex hull", y = "Cumulative species richness") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

# Save session
save.image("output/sSBR_comparison.RData")
