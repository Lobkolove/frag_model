library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

source("R/prepare_data_scale_new.R")
source("R/prepare_data_disp_new.R")
source("R/prepare_data_heatmap_disp.R")
source("R/prepare_data_heatmap_disp_sample.R")

palette_dis4 <- c("#D50A0A", "#013369", "#f97676", "#3796fd", "#006F3D", "#26FD9B")

# Figure 3 ------------------------------------------------------------------

cur_path <- "data-raw/model_output/fig_3"

cur_sim <- prepare_data_scale_new(cur_path)

landscape <- cur_sim$landscape
sample <- cur_sim$sample

panel1 <- ggplot(data = landscape, aes(x = as.factor(fragmentation), y = richness_val, color = as.factor(richness_scale))) +
  stat_summary(fun.data = "mean_se", geom = "point", size = 3) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, aes(color = as.factor(richness_scale)), size = 1.2) +
  labs(x = "fragmentation level", y = "species richness", title = "Landscape scale") +
  scale_color_discrete(type = palette_dis4[c(1, 3)], name = "", labels = c("Net effects", "Geometric effects")) +
  geom_smooth(method = lm, se = FALSE, aes(group = as.factor(richness_scale), color = as.factor(richness_scale)), size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "bottom")

panel2 <- ggplot(data = sample, aes(x = as.factor(fragmentation), y = richness_val, color = as.factor(richness_scale))) +
  stat_summary(fun.data = "mean_se", geom = "point", size = 3) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, aes(color = as.factor(richness_scale)), size = 1.2) +
  labs(x = "fragmentation level", y = "species richness", title = "Sample scale") +
  scale_color_discrete(type = palette_dis4[c(2, 4)], name = "", labels = c("Net effects", "Geometric effects")) +
  geom_smooth(method = lm, se = FALSE, aes(group = as.factor(richness_scale), color = as.factor(richness_scale)), size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "bottom")

# Function for extracting legend only
get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

# Extract legends
legend1 <- get_only_legend(panel1)
legend2 <- get_only_legend(panel2)

# Remove legends from panel1 and panel2
panel1 <- panel1 + theme(legend.position = "none")
panel2 <- panel2 + theme(legend.position = "none")

# Add labels to panels using cowplot
panel1 <- ggdraw(panel1) + draw_label("(b)", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 15)
panel2 <- ggdraw(panel2) + draw_label("(a)", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 15)

# Arrange the plots
lay <- rbind(
  c(1, 2),
  c(4, 3)
)

hig <- c(9, 1)

multi_plot <- gridExtra::grid.arrange(panel2, panel1, legend1, legend2, nrow = 2, ncol = 3, heights = hig, layout_matrix = lay)

# Save the figure
ggsave("R/figures/fig3.jpg", plot = multi_plot, width = 12, height = 6, dpi = 400)

# Figure 4 ------------------------------------------------------------------

cur_path <- "data-raw/model_output/fig_4"

current_sim <- prepare_data_disp_new(cur_path)
cur_sim <- current_sim$data

disp_levels <- c(1, 3, 5)


data_disp1 <- cur_sim %>%
  filter(disp_dist == disp_levels[1])

plot_lan_1 <- ggplot(data_disp1, aes(x = as.factor(fragmentation), y = present_species, color = as.factor(step))) +
  stat_summary(fun.data = "mean_se", geom = "point", size = 1.5) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, aes(color = as.factor(step)), size = 0.6) +
  theme_bw() +
  scale_color_discrete(type = palette_dis4[c(3, 1)], name = "", labels = c("Geometric", "Net effects")) +
  labs(x = "Fragmentation", y = "Richness") +
  geom_smooth(method = "lm", se = FALSE, aes(group = step), linewidth = 0.5) +
  theme(
    legend.position = "none",
    axis.title = element_blank()
  ) +
  ylim(0, 900)

data_disp2 <- cur_sim %>%
  filter(disp_dist == disp_levels[2])

plot_lan_2 <- ggplot(data_disp2, aes(x = as.factor(fragmentation), y = present_species, color = as.factor(step))) +
  stat_summary(fun.data = "mean_se", geom = "point", size = 1.5) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, aes(color = as.factor(step)), size = 0.6) +
  theme_bw() +
  scale_color_discrete(type = palette_dis4[c(3, 1)], name = "", labels = c("Geometric", "Net effects")) +
  labs(x = "Fragmentation", y = "Richness") +
  geom_smooth(method = "lm", se = FALSE, aes(group = step), linewidth = 0.5) +
  theme(
    legend.position = "none",
    axis.title = element_blank()
  ) +
  ylim(0, 900)

data_disp3 <- cur_sim %>%
  filter(disp_dist == disp_levels[3])

plot_lan_3 <- ggplot(data_disp3, aes(x = as.factor(fragmentation), y = present_species, color = as.factor(step))) +
  stat_summary(fun.data = "mean_se", geom = "point", size = 1.5) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, aes(color = as.factor(step)), size = 0.6) +
  theme_bw() +
  scale_color_discrete(type = palette_dis4[c(3, 1)], name = "", labels = c("Geometric", "Net effects")) +
  labs(x = "Fragmentation", y = "Richness") +
  geom_smooth(method = "lm", se = FALSE, aes(group = step), linewidth = 0.5) +
  theme(axis.title = element_blank()) +
  ylim(0, 900)

plot_sam_1 <- ggplot(data_disp1, aes(x = as.factor(fragmentation), y = sample_species_mean, color = as.factor(step))) +
  stat_summary(fun.data = "mean_se", geom = "point", size = 1.5) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, aes(color = as.factor(step)), size = 0.6) +
  theme_bw() +
  scale_color_discrete(type = palette_dis4[c(4, 2)], name = "", labels = c("Geometric", "Net effects")) +
  labs(x = "Fragmentation", y = "Richness") +
  geom_smooth(method = "lm", se = FALSE, aes(group = step), linewidth = 0.5) +
  theme(
    legend.position = "none",
    axis.title = element_blank()
  ) +
  ylim(0, 30)

plot_sam_2 <- ggplot(data_disp2, aes(x = as.factor(fragmentation), y = sample_species_mean, color = as.factor(step))) +
  stat_summary(fun.data = "mean_se", geom = "point", size = 1.5) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, aes(color = as.factor(step)), size = 0.6) +
  theme_bw() +
  scale_color_discrete(type = palette_dis4[c(4, 2)], name = "", labels = c("Geometric", "Net effects")) +
  labs(x = "Fragmentation", y = "Richness") +
  geom_smooth(method = "lm", se = FALSE, aes(group = step), linewidth = 0.5) +
  theme(
    legend.position = "none",
    axis.title = element_blank()
  ) +
  ylim(0, 30)

plot_sam_3 <- ggplot(data_disp3, aes(x = as.factor(fragmentation), y = sample_species_mean, color = as.factor(step))) +
  stat_summary(fun.data = "mean_se", geom = "point", size = 1.5) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2, aes(color = as.factor(step)), size = 0.6) +
  theme_bw() +
  scale_color_discrete(type = palette_dis4[c(4, 2)], name = "", labels = c("Geometric", "Net effects")) +
  labs(x = "Fragmentation", y = "Richness") +
  geom_smooth(method = "lm", se = FALSE, aes(group = step), linewidth = 0.5) +
  theme(axis.title = element_blank()) +
  ylim(0, 30)


# aggregated figure

data_lm <- current_sim$data_lm

# add demo value at each dispersal distance to data frame

disp_vals <- unique(data_lm$disp_dist)
steps <- unique(data_lm$step)
dem_val <- list()
for (i in 1:6) {
  obj <- data_lm$estimate[data_lm$step == steps[1] & data_lm$disp_dist == disp_vals[i]] - data_lm$estimate[data_lm$step == steps[2] & data_lm$disp_dist == disp_vals[i]]
  dem_val <- append(dem_val, obj)
}


agg <- ggplot(data_lm, aes(x = disp_dist, y = estimate, color = as.factor(step))) +
  geom_point(size = 5) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  labs(x = "Dispersal distance", y = "Slope", title = "Aggregated") +
  scale_color_discrete(type = palette_dis4[c(3, 1)], name = "", labels = c("Geometric", "Net effects")) +
  scale_x_continuous("Dispersal distance", labels = c(1:6), breaks = c(1:6)) +
  geom_hline(yintercept = 0, linetype = "dotted", col = "black") +
  geom_segment(
    x = 1, y = 515,
    xend = 1, yend = 290,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 1,
    arrow = arrow(length = unit(0.1, "inches"), ends = "both"),
    colour = "darkgreen"
  ) +
  geom_segment(
    x = 2, y = 450,
    xend = 2, yend = 115,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 1,
    arrow = arrow(length = unit(0.1, "inches"), ends = "both"),
    colour = "darkgreen"
  ) +
  geom_segment(
    x = 3, y = 400,
    xend = 3, yend = -10,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 1,
    arrow = arrow(length = unit(0.1, "inches"), ends = "both"),
    colour = "darkgreen"
  ) +
  geom_segment(
    x = 4, y = 400,
    xend = 4, yend = -100,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 1,
    arrow = arrow(length = unit(0.1, "inches"), ends = "both"),
    colour = "darkgreen"
  ) +
  geom_segment(
    x = 5, y = 330,
    xend = 5, yend = -155,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 1,
    arrow = arrow(length = unit(0.1, "inches"), ends = "both"),
    colour = "darkgreen"
  ) +
  geom_segment(
    x = 6, y = 270,
    xend = 6, yend = -155,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 1,
    arrow = arrow(length = unit(0.1, "inches"), ends = "both"),
    colour = "darkgreen"
  ) +
  annotate("segment",
    x = 6.5, y = 20, xend = 6.5, yend = -40,
    col = "darkgreen", arrow = arrow(length = unit(0.08, "inches"), ends = "both")
  ) +
  annotate("text",
    x = 6.91, y = -10,
    label = "Demographic", size = 3
  ) +
  annotate("text",
    x = 1.1, y = round(dem_val[[1]]) / 2 + data_lm$estimate[data_lm$step == steps[2] & data_lm$disp_dist == 1],
    label = paste(round(dem_val[[1]]) * -1), size = 3
  ) +
  annotate("text",
    x = 2.1, y = round(dem_val[[2]]) / 2 + data_lm$estimate[data_lm$step == steps[2] & data_lm$disp_dist == 2],
    label = paste(round(dem_val[[2]]) * -1), size = 3
  ) +
  annotate("text",
    x = 3.1, y = round(dem_val[[3]]) / 2 + data_lm$estimate[data_lm$step == steps[2] & data_lm$disp_dist == 3],
    label = paste(round(dem_val[[3]]) * -1), size = 3
  ) +
  annotate("text",
    x = 4.1, y = round(dem_val[[4]]) / 2 + data_lm$estimate[data_lm$step == steps[2] & data_lm$disp_dist == 4],
    label = paste(round(dem_val[[4]]) * -1), size = 3
  ) +
  annotate("text",
    x = 5.1, y = round(dem_val[[5]]) / 2 + data_lm$estimate[data_lm$step == steps[2] & data_lm$disp_dist == 5],
    label = paste(round(dem_val[[5]]) * -1), size = 3
  ) +
  annotate("text",
    x = 6.1, y = round(dem_val[[6]]) / 2 + data_lm$estimate[data_lm$step == steps[2] & data_lm$disp_dist == 6],
    label = paste(round(dem_val[[6]]) * -1), size = 3
  ) +
  coord_cartesian(xlim = c(1, 6), clip = "off")

tit_1 <- textGrob("Dispersal distance = 1", gp = gpar(fontsize = 14))
tit_2 <- textGrob("Dispersal distance = 3", gp = gpar(fontsize = 14))
tit_3 <- textGrob("Dispersal distance = 5", gp = gpar(fontsize = 14))
tit_lan <- textGrob("Landscape", gp = gpar(fontsize = 13), rot = 90)
tit_sam <- textGrob("Sample", gp = gpar(fontsize = 13), rot = 90)
tit_spec <- textGrob("Species richness", gp = gpar(fontsize = 15), rot = 90)
tit_frag <- textGrob("Fragmentation level", gp = gpar(fontsize = 15))

# Add labels to each plot
plot_lan_1 <- plot_lan_1 + ggtitle("(d)")
plot_lan_2 <- plot_lan_2 + ggtitle("(e)")
plot_lan_3 <- plot_lan_3 + ggtitle("(f)")
plot_sam_1 <- plot_sam_1 + ggtitle("(a)")
plot_sam_2 <- plot_sam_2 + ggtitle("(b)")
plot_sam_3 <- plot_sam_3 + ggtitle("(c)")
agg <- agg + ggtitle("(g)")

# Adjust the layout to include the titles
layout <- rbind(
  c(NA, NA, 1, 2, 3, NA),
  c(13, 4, 5, 6, 7, 7),
  c(13, 8, 9, 10, 11, 11),
  c(NA, NA, NA, 14, NA, NA),
  c(NA, 12, 12, 12, 12, 12)
)

# Create the multi-panel plot with labels
mp <- grid.arrange(
  tit_1, tit_2, tit_3, tit_sam, plot_sam_1, plot_sam_2, plot_sam_3,
  tit_lan, plot_lan_1, plot_lan_2, plot_lan_3, agg, tit_spec, tit_frag,
  layout_matrix = layout,
  widths = c(.05, .05, .265, .265, .27, .1),
  heights = c(.05, .25, .25, .05, .4)
)

ggsave("R/figures/fig4.jpg", plot = mp, width = 10, height = 13.5, dpi = 400)

# Figure 5 ------------------------------------------------------------------

cur_path <- "data-raw/model_output/fig_5"

cur_sim <- prepare_data_heatmap_disp(cur_path, "disp_dist", "edge")

geo <- cur_sim$geo_lm
gedo <- cur_sim$gedo_lm
demo <- as.data.frame(geo) %>%
  select(-c(data, statistic, p.value, std.error)) %>%
  mutate(estimate = gedo$estimate - geo$estimate)

# get total simulations and extinctions

extinct <- cur_sim$data %>%
  filter(step == cur_sim$steps[2]) %>%
  tally(individuals == 0) %>%
  as.numeric()

tot_sim <- cur_sim$data %>%
  filter(step == cur_sim$steps[2]) %>%
  tally() %>%
  as.numeric()

p1 <- ggplot(cur_sim$gedo_lm, aes(as.factor(disp_dist), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Net fragmentation effects", x = "dispersal distance", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "near", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "far", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(c)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

p2 <- ggplot(cur_sim$geo_lm, aes(as.factor(disp_dist), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Geometric fragmentation effects", x = "dispersal distance", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "near", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "far", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(a)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

p3 <- ggplot(demo, aes(as.factor(disp_dist), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Demographic fragmentation effects", x = "dispersal distance", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "near", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "far", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(b)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

multi_plot <- gridExtra::grid.arrange(p2, p3, p1, ncol = 1)
ggsave("R/figures/fig5.jpg", plot = multi_plot, width = 6, height = 13, dpi = 400)

# Figure 6 ------------------------------------------------------------------

cur_path <- "data-raw/model_output/fig_6"

cur_sim <- prepare_data_heatmap_disp(cur_path, "ac", "edge")

extinct <- cur_sim$data %>%
  filter(step == cur_sim$steps[2]) %>%
  tally(individuals == 0) %>%
  as.numeric()

tot_sim <- cur_sim$data %>%
  filter(step == cur_sim$steps[2]) %>%
  tally() %>%
  as.numeric()

geo <- cur_sim$geo_lm
gedo <- cur_sim$gedo_lm
demo <- as.data.frame(geo) %>%
  select(-c(data, statistic, p.value, std.error)) %>%
  mutate(estimate = gedo$estimate - geo$estimate)

p3 <- ggplot(demo, aes(as.factor(ac), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Demographic fragmentation effects", x = "landscape autocorrelation", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "low", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "high", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(b)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")


p2 <- ggplot(cur_sim$geo_lm, aes(as.factor(ac), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Geometric fragmentation effects", x = "landscape autocorrelation", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "low", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "high", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(a)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

p1 <- ggplot(cur_sim$gedo_lm, aes(as.factor(ac), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Net fragmentation effects", x = "landscape autocorrelation", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "low", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "high", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(c)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

multi_plot <- gridExtra::grid.arrange(p2, p3, p1, ncol = 1)
ggsave("R/figures/fig6.jpg", plot = multi_plot, width = 6, height = 13, dpi = 400)



# Figure S3A ----------------------------------------------------------------

cur_path <- "data-raw/model_output/fig_s3a"

cur_sim <- prepare_data_heatmap_disp_samp(cur_path, "disp_dist", "edge")

geo <- cur_sim$geo_lm_sam
gedo <- cur_sim$gedo_lm_sam
demo <- as.data.frame(geo) %>%
  select(-c(data, statistic, p.value, std.error)) %>%
  mutate(estimate = gedo$estimate - geo$estimate)

p1 <- ggplot(cur_sim$gedo_lm_sam, aes(as.factor(disp_dist), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Net effects", x = "dispersal distance", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "near", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "far", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(c)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

p2 <- ggplot(cur_sim$geo_lm_sam, aes(as.factor(disp_dist), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Geometric", x = "dispersal distance", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "near", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "far", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(c)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

p3 <- ggplot(demo, aes(as.factor(disp_dist), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Demographic", x = "dispersal distance", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "near", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "far", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(a)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

multi_plot <- gridExtra::grid.arrange(p2, p3, p1, ncol = 1)
ggsave("R/figures/figs3a.png", plot = multi_plot, width = 6, height = 13)

# Figure S3B ----------------------------------------------------------------

cur_path <- "data-raw/model_output/fig_s3b"

cur_sim <- prepare_data_heatmap_disp_samp(cur_path, "ac", "edge")

geo <- cur_sim$geo_lm_sam
gedo <- cur_sim$gedo_lm_sam
demo <- as.data.frame(geo) %>%
  select(-c(data, statistic, p.value, std.error)) %>%
  mutate(estimate = gedo$estimate - geo$estimate)
# get total simulations and extinctions

extinct <- cur_sim$data %>%
  filter(step == cur_sim$steps[2]) %>%
  tally(individuals == 0) %>%
  as.numeric()

tot_sim <- cur_sim$data %>%
  filter(step == cur_sim$steps[2]) %>%
  tally() %>%
  as.numeric()

p1 <- ggplot(cur_sim$gedo_lm_sam, aes(as.factor(ac), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Net effects", x = "landscape autocorrelation", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "near", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "far", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(c)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

p2 <- ggplot(cur_sim$geo_lm_sam, aes(as.factor(ac), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Geometric", x = "landscape autocorrelation", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "near", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "far", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(a)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

p3 <- ggplot(demo, aes(as.factor(ac), as.factor(edge), fill = estimate)) +
  geom_tile() +
  scale_fill_gradient2(low = "coral2", mid = "white", high = "steelblue", midpoint = 0) +
  labs(title = "Demographic", x = "landscape autocorrelation", y = "edge effects", fill = "Slope") +
  # geom_label(aes(label = round(estimate, 1)), fill = "white" , size = 2.5) +
  theme(panel.background = element_blank(), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  annotate("text", x = 1, y = -0.6, label = "near", size = 3.5) +
  annotate("text", x = 10, y = -0.6, label = "far", size = 3.5) +
  annotate("text", x = -0.5, y = 1.6, label = "positive", size = 3.5, angle = 90) +
  annotate("text", x = -0.5, y = 9.3, label = "negative", size = 3.5, angle = 90) +
  annotate("text", x = 11.5, y = 10, label = "(b)", size = 8) +
  coord_cartesian(ylim = c(1, 10), xlim = c(1, 10), clip = "off")

multi_plot <- gridExtra::grid.arrange(p2, p3, p1, ncol = 1)
ggsave("R/figures/figs3b.png", plot = multi_plot, width = 6, height = 13)
