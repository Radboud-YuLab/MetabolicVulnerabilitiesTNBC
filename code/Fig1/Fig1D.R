## 3D growht curves plto
# written by: Cyriel Huijer 
# Last updated: 22 June 2026

# Load packages
library(tidyverse)

# Clear environment
rm(list=ls())

# Set working directory to /code
setwd("")


# Load data
data <- read.csv("data/growth_rate/celltiterglo_3D/cell_count_3D_feed_R.csv")
data$condition_linked <- paste0(data$condition, "_", data$conc)

# THe order of bars
order_vec <- c(
  "EV_0", "EV_1000",
  "MYC_0","MYC_80","MYC_160",  "MYC_1000",
  paste0("TP53R175H_MYC_", c(0,80,160,1000)),
  paste0("TP53R248Q_MYC_", c(0,80,160,1000)),
  paste0("TP53R273H_MYC_", c(0,80,160,1000))
)

data

# Calculate growth rate
slopes <- data %>%
  group_by(condition_linked, rep) %>%
  do({
    fit <- lm(ln_cells ~ timepoint, data = .)
    # extract slope
    slope <- coef(fit)[["timepoint"]]
    # extract R2
    r2 <- summary(fit)$r.squared
    data.frame(slope = slope, r_squared = r2)
  }) %>%
  ungroup()

slopes <- slopes %>%
  filter(r_squared > 0.4) %>%
  filter(rep != 2)

slopes_avg <- data %>%
  group_by(condition_linked) %>%
  do({
    fit <- lm(ln_cells ~ timepoint, data = .)
    # extract slope
    slope <- coef(fit)[["timepoint"]]
    # extract R2
    r2 <- summary(fit)$r.squared
    data.frame(slope = slope, r_squared = r2)
  }) %>%
  ungroup()



# Cells appear to be in exponential growth phase between 24 and 72 hours.
# Not enough datapoints to accurately calculate slope and R2. GR is based on t24 and t72

gr_data <- data %>%
  filter(timepoint %in% c(24,72))


growth_rates <- gr_data %>%
  group_by(condition_linked, rep) %>%
  summarise(
    n24 = ln_cells[timepoint == 24],
    n72 = ln_cells[timepoint == 72],
    k = (n72 - n24) / (72 - 24),
    .groups = "drop"
  )

growth_rates <- growth_rates %>%
  filter(k > 0 )

growth_rates <- growth_rates %>%
  mutate(subline = str_extract(condition_linked, ".*(?=_[^_]+$)"))

# File for exch. flux calculation:
#write.csv(growth_rates, "data/growth_rate/celltiterglo_3D/growth_rates_filtered_reps_no_filter.csv")

growth_rates_filtered <- growth_rates

#  color palette
colors <- c(
  "TP53R175H_MYC" = "#D81B60", 
  "TP53R248Q_MYC" = "#1E88E5", 
  "TP53R273H_MYC" = "#FFC107",
  "MYC" = "#004D40", 
  "EV" = "#4EC6A9"
)

growth_rates_filtered$condition_linked <- factor(
  growth_rates_filtered$condition_linked,
  levels = order_vec
)

#write.csv(growth_rates_filtered, "../../../1_manuscript/files_figures/growth_rates_filtered_reps.csv")

# calculate the means of growth rate
means_k <- growth_rates_filtered %>%
  group_by(condition_linked, subline) %>%
  summarise(mean_k = mean(k), .groups = "drop") %>%
  filter(!condition_linked %in% c("MYC_80","MYC_160"))
growth_rates_filtered <- growth_rates_filtered %>%
  filter(!condition_linked %in% c("MYC_80","MYC_160"))


# Plot Fig 1D
fig1D <- ggplot(growth_rates_filtered,
                      aes(x = condition_linked, y = k)) +
  
  # barplot of mean
  geom_col(data = means_k,
           aes(x = condition_linked, y = mean_k, fill = subline),
           width = 0.8,
           alpha = 1) +
  
  # replicates
  geom_point(size = 3, alpha = 0.9) +
  
  scale_fill_manual(values = colors) +
  
  labs(
       y = "Growth rate",
       x = "Condition") +
  
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
  
fig1D

#ggsave("../../../1_manuscript/figures/final_figures/Fig3/growth_rate_3D_bar_no_filter.pdf", plot = plot_growth, width = (16/12)*5, height = 3, units = "in", dpi = 300)
