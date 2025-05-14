# create figure 4
# show the SAR for 10 randomly selected plots

full_dob_neon_richness_data <- readRDS("~/soil-sar/plot-sar-permutation/full_dob_neon_richness_data.rds")
load("~/soil-sar/plot-sar-permutation/standard sar neon/richness_subplot10_neon_standar.RData")
full_dob_neon_richness_data_rarefaction <- full_dob_neon_richness_data

full_dob_neon_richness_data_rarefaction %>%
  filter(area != 100 & guild == "all" & !is.na(sd_value)) %>%
  dplyr::select(plotid, mean_value, sd_value, area) %>%
  bind_rows(richness_subplot10_neon_standar) -> full_neon_data_with_sd


full_neon_data_with_sd %>%
  group_by(plotid) %>%
  summarize(n()) %>%
  filter(`n()` > 3) %>%
  pull(plotid) -> plotid_with_four_point

set.seed(125)
rand_plot <- sample(plotid_with_four_point, 10, replace = FALSE)
full_neon_data_with_sd %>% filter(plotid %in% rand_plot & !is.na(mean_value)) -> plot_data



# the plot with the raw data
# Fit a nonlinear least squares model for each group
fit_power_law <- function(df) {
  nls(mean_value ~ c * area^z, data = df, start = list(c = 1, z = 0.5))
}

# Apply the NLS function to each group
nls_fits <- plot_data %>%
  group_by(plotid) %>%
  do(model = fit_power_law(.))

# to get the values

nls_fits$model[[1]]

# Predict values based on the model fits
predicted_values <- nls_fits %>%
  rowwise() %>%
  mutate(Predicted = list(data.frame(
    area = seq(min(plot_data$area), max(plot_data$area), length.out = 100),
    mean_value = predict(model, newdata = data.frame(area = seq(min(plot_data$area), max(plot_data$area), length.out = 100)))
  ))) %>%
  unnest(cols = c(Predicted))

# Plotting the data and the fitted NLS lines

p_zvalue <- ggplot(plot_data, aes(x = area, y = mean_value, color = plotid)) +
  geom_errorbar(data = plot_data, aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, color = plotid), width = 45) +
  geom_point(size = 3) +
  scale_color_manual("",
    breaks = unique(plot_data$plotid),
    labels = c(
      expression("Central Plains Experimental Range " * italic(z) * "=0.65"),
      expression("Guanica Dry Forest Reserve " * italic(z) * "=0.74"),
      expression("Harvard Forest " * italic(z) * "=0.65"),
      expression("Niwot Ridge Mountain Research Station " * italic(z) * "=0.70"),
      expression("Oak Ridge National Laboratory " * italic(z) * "=0.69"),
      expression("Ordway-Swisher Biological Station " * italic(z) * "=0.82"),
      expression("The San Joaquin Experimental range " * italic(z) * "=0.65"),
      expression("The Talladega National Forest " * italic(z) * "=0.79"),
      expression("University of Kansas Field Station " * italic(z) * "=0.72"),
      expression("Woodworth " * italic(z) * "=0.66")
    ),
    values = c("chocolate1", "#037f77", "royalblue", "#f0a73a", "forestgreen", "#7c1a97", "#c94e65", "tan", "pink", "gray")
  ) +
  geom_line(data = predicted_values, aes(x = area, y = mean_value, color = plotid), linetype = "solid", size = 1.2) +
  theme(
    legend.position = c(0.37, 0.780361205205),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.key.height = unit(0.48, "cm"),
    text = element_text(size = 18),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_text(hjust = 0, size = 15),
    axis.text.x = element_text(hjust = 0.5, size = 15),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    plot.margin = unit(c(0.2, 2, 0, 1), "cm"),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1, fill = NA)
  ) +
  xlab(expression(Area * " " * (m^2))) +
  ylab("Taxonomic diversity") +
  guides(color = guide_legend(nrow = 10)) +
  ylim(0, 3500)




p_z_guild <- ggplot(guild_mean, aes(x = guild, y = zvalue, fill = guild), alpha = 0.5) +
  # sm_raincloud(size=0.1,point.params =list(size=2),sep_level=2)+
  geom_boxplot(data = guild_mean, aes(x = guild, y = zvalue, fill = guild), width = 0.4, size = 0.5, color = "black", size = 0.1, outlier.size = 1) +
  scale_fill_manual("",
    breaks = od$guild, values = c("chocolate1", "#037f77", "royalblue", "#c94e65"),
    labels = c(
      "AM (N=326)", "EM (N=415)",
      "Plant pathogens (N=403)",
      "Soil saprotrophs (N=416)"
    )
  ) +
  theme(
    legend.position = c(0.5, 0.12660),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    text = element_text(size = 18),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.ticks.x = element_blank(),
    # axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  ylab(expression(italic(z) * " value")) +
  xlab("Trophic guilds") +
  geom_hline(yintercept = 0.71, color = "red", linetype = "dashed", size = 0.8) +
  annotate("text", x = 1, y = 1.3, label = "a", size = 6) +
  annotate("text", x = 2, y = 1, label = "a", size = 6) +
  annotate("text", x = 3, y = 1.2, label = "b", size = 6) +
  annotate("text", x = 4, y = 1, label = "c", size = 6) +
  # annotate("text", x = 5, y = 1.2, label = "d", size = 6) +
  # annotate("text", x = 6, y = 0.95, label = "e", size = 6) +
  # annotate("text", x = 7, y = 1.014, label = "d", size = 6) +
  # annotate("text", x = 8, y = 1.15, label = "d", size = 6) +
  ylim(0, 1.3)



# display the distribution of the estimated z values with a histgram

p_z_distribu <- ggplot(full_parameter_data %>% filter(guild == "all" & !is.na(zvalue)), aes(x = zvalue)) +
  geom_histogram(binwidth = 0.02, fill = "#F3A332", color = "black", alpha = 0.7) +
  ggtitle("") +
  xlab(expression(italic(z) * " value")) +
  ylab("Frequency") +
  geom_vline(xintercept = 0.713, linetype = "dashed", color = "red") +
  theme(
    legend.position = c(0.75, 0.28),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    text = element_text(size = 18),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_text(hjust = 0, size = 15),
    axis.text.x = element_text(hjust = 0.5, size = 15),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    plot.margin = unit(c(0.2, 0.5, 0, 0.5), "cm"),
    panel.background = element_rect(fill = "NA"),
    panel.border = element_rect(color = "black", size = 1, fill = NA)
  ) +
  annotate("text", x = 0.5, y = 50, label = "CV=11%", size = 5)


p11 <- ggplot(data = effect_no_plant %>% filter(guild == "all"), aes(x = estimate, y = 1:13)) +
  geom_point(size = 4, pch = 21, color = "black", fill = rep(c("#DE582E", "#F3A332", "#1868b2"), times = c(1, 5, 7))) +
  geom_errorbar(aes(xmin = estimate - 1.96 * sd, xmax = 1.96 * sd + estimate, width = 0.2),
    color = rep(c("#DE582E", "#F3A332", "#1868b2"), times = c(1, 5, 7))
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_continuous(breaks = 1:13, labels = c(
    expression(italic(C)), "pH", "Moisture", "SoilC", "CEC", "Sand", "MAT",
    "MDR", "Temp.seas.", "MTWQ", "MAP", "Pre.seas.", "PWQ"
  )) +
  theme(
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 15),
    axis.ticks.length.y = unit(-0.150, "cm"),
    axis.text.y = element_text(size = 12, hjust = 0, margin = margin(r = -70)),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.x = element_text(angle = 0, size = 15, hjust = 0.5),
    panel.background = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.3), "cm")
  ) +
  ylab("") +
  xlab("Effect size ± 95% CI (N=453)") +
  xlim(-0.6, 0.6) +
  # ggtitle("Climate+Soil(N=453)")+
  annotate("text", x = -0.31, y = 1, label = "***", size = 8) +
  annotate("text", x = 0.265, y = 7, label = "*", size = 8) +
  annotate("text", x = 0.265, y = 9, label = "*", size = 8) +
  annotate("text", x = 0.24583, y = 11, label = "*", size = 8) +
  annotate("text", x = -0.16820583, y = 3, label = "*", size = 8) +
  xlim(-0.5, 0.5)


p1 <- ggplot(data = varp_new_noplant %>% filter(Fractions > 0 & guild == "all"), aes(x = guild, y = Fractions, fill = type)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.6) +
  scale_fill_manual(expression(italic(R["Fixed"]^2) * " =58.14%"),
    breaks = unique(varp_new_noplant$type),
    labels = c(
      expression(italic(C)), "S", "Cli.", expression(italic(C) * "+S"), expression(italic(C) * "+Cli."),
      "S+Cli.", expression(italic(C) * "+S+Cli."), "Random", "Residuals"
    ),
    values = c("#F3A332", "seagreen1", "royalblue", "greenyellow", "forestgreen", "purple", "lavender", "#FFFFCC", "gray")
  ) +
  theme(
    legend.position = c(0.5, 0.7),
    panel.border = element_rect(fill = NA, size = 1, color = "black"),
    axis.title.x = element_text(size = 15),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 18),
    legend.key.size = unit(0.35, "cm"),
    legend.key.height = unit(0.38, "cm"),
    axis.text.y = element_text(size = 15, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm"),
    panel.background = element_rect(fill = "NA")
  ) +
  guides(fill = guide_legend(nrow = 10, byrow = TRUE)) +
  ylab("Variance explained") +
  xlab("") +
  # ggtitle(expression(italic(R["Fixed"]^2)*" =56.4%"))+
  ylim(0, 2)



p1 <- ggplotGrob(p1)
p11 <- ggplotGrob(p11)

p1$heights <- p11$heights

p_effects <- plot_grid(p1, p11, ncol = 2, rel_widths = c(1, 2))


p_z_pva <- ggplotGrob(p_zvalue)
p_effects <- ggplotGrob(p_effects)
p_z_guild <- ggplotGrob(p_z_guild)
p_z_distribu <- ggplotGrob(p_z_distribu)
# align the heights
p_z_pva$heights <- p_z_distribu$heights
p_z_pva$heights <- p_z_guild$heights
p_z_pva$heights <- p_z_distribu$heights
p_z_guild$heights <- p_effects$heights


p_z_pva$widths <- p_z_distribu$widths
p_z_pva$widths <- p_z_guild$widths
p_z_guild$widths <- p_effects$widths
p_z_guild$widths <- p_z_pva$widths

p_z_guild$heights <- p_z_distribu$heights # important to align the heights

plot_grid(p_z_pva, p_z_distribu, p_z_guild, p_effects,
  ncol = 2, label_size = 18, label_x = 0.1,
  label_y = c(1.01, 1.01, 1.03, 1.03),
  labels = paste0("(", letters[1:4], ")")
)

plot_grid(p_z_guild, p_z_distribu)
