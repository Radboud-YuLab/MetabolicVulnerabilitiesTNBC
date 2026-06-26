# Growth Curves Double Oncogene Lines (2D):
# Programmer: Cyriel Huijer
# Date: 22 May 2024
# Updated: 22 June 2022

# Load in libraries:
library(tidyverse)
library(reshape2)
library(dplyr)
library(pracma)
library(broom)

# Clear environment
rm(list=ls())

# Set working directory to /code
setwd("")

#Load files 
rep1 <- read.csv("data/growth_rate/IncuCyte_2D/CH00022_doubleOncogenes_rep1.csv")
rep2 <- read.csv("data/growth_rate/IncuCyte_2D/CH00024_doubleOncogenes_rep2.csv")

# format data
melted_df_1 <- melt(rep1, id.vars = "timepoint", variable.name = "condition", value.name = "value")
melted_df_1 <- melted_df_1 %>%
  separate(condition, into = c("subline", "dox_concentration", "replicate"), sep = "_")
mean_values_1 <- melted_df_1 %>%
  group_by(subline, dox_concentration,timepoint) %>%
  summarise(mean_value = mean(value))
mean_values_1$exp <- 1


melted_df_2 <- melt(rep2, id.vars = "timepoint", variable.name = "condition", value.name = "value")
melted_df_2 <- melted_df_2 %>%
  separate(condition, into = c("subline", "dox_concentration", "replicate"), sep = "_")
melted_df_2 <- subset(melted_df_2,timepoint <= 48)
mean_values_2 <- melted_df_2 %>%
  group_by(subline, dox_concentration,timepoint) %>%
  summarise(mean_value = mean(value))
mean_values_2$exp <- 2

combined <- rbind(mean_values_1, mean_values_2)
combined$subline <- factor(combined$subline, levels = c("EV", "MYC", "TP53R175H", "TM"))

# Cell count comes from the following file which calculates the cell count/confluency:IncuCyte_CyQUANT_benchmark.xlsx

# conversion of confluency to cell count
combined$cell_count <- (3623.7*exp(0.0271*combined$mean_value))#*(9.6/0.32)
combined$log_cell_count <- log(combined$cell_count)
# extract datapoints in log phase
log_phase <- combined %>%
  subset(timepoint > 13) %>%
  subset(timepoint < 35)

# calculate growth rate
slopes <- log_phase %>%
  group_by(subline,dox_concentration,exp) %>%
  do(tidy(lm(log_cell_count ~ timepoint, data = .))) %>%
  filter(term == "timepoint") %>%
  dplyr::select(subline, slope = estimate)
mean(slopes$slope)
  
# Colors of cell lines
colors <- c(
  "MYC" = "#004D40",
  "EV"  = "#4EC6A9",  
  "TP53R175H" = "#6A0E97",
  "TM" = "#D81B60"  
)

# THe order of bars
order_vec <- c(
  "EV_0", "EV_100","EV_1000",
  "MYC_0", "MYC_100","MYC_1000",
  "TP53R175H_0","TP53R175H_100","TP53R175H_1000",
  "TM_0","TM_100","TM_1000"
)


means_slopes <- slopes %>%
  group_by(subline, dox_concentration) %>%
  summarize(mean_slope = mean(slope), .groups = "drop")

# create combined label for ordering
slopes$sub_dox <- paste0(slopes$subline, "_", slopes$dox_concentration)
means_slopes$sub_dox <- paste0(means_slopes$subline, "_", means_slopes$dox_concentration)

# ordering with factors
slopes$sub_dox <- factor(slopes$sub_dox, levels = order_vec)
means_slopes$sub_dox <- factor(means_slopes$sub_dox, levels = order_vec)

p1 <- ggplot() +
  geom_col(data = means_slopes,
           aes(x = sub_dox,
               y = mean_slope,
               fill = subline),
           width = 0.8) +
  geom_point(data = slopes,
             aes(x = sub_dox,
                 y = slope),
             size = 2,
             fill = "black",
             color = "black",
             stroke = 0.8) +
  # add colors
  scale_fill_manual(values = colors) +
  # hide legend
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "Growth Rate (h-1)", x = "")
p1

# Plot fig 1C
#ggsave("../../../1_manuscript/figures/final_figures/Fig3/growth_rate_2D_bar.pdf", plot = p1, width = 5, height = 3, units = "in", dpi = 300)
