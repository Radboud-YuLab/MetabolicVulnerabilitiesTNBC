## Percent changed subsystem plots (Fig 3B):
# written by: Cyriel Huijer
# last updated: 22 June 2026

library(tidyverse)
library(matrixStats)

# clear environment
rm(list=ls())

# Set working directory to /code
setwd("")

# Load in dataframe from MATLAB

data <- read.csv("data/modeling/output/20251218_pct_changed_subsystems.csv")
data$median <- matrixStats::rowMedians(as.matrix(data[, 2:6]))

data <- data[order(data$median,decreasing = TRUE),]
data<- data[1:10,]

data_long <- pivot_longer(data, cols = 2:6, names_to = "Condition", values_to = "Pct_Change")
data_long$subsytem <- factor(data_long$subsytem, levels = rev(data$subsytem))

colors <- c(
  "TP53R175H_MYC" = "#D81B60", 
  "TP53R248Q_MYC" = "#1E88E5",
  "TP53R273H_MYC" = "#FFC107",
  "MYC" = "#004D40", 
  "EV" = "#4EC6A9"
)

# Fig 3B
ggplot(data_long, aes(x = Pct_Change, y = subsytem, fill = Condition)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs(title = "% rxn change per subsystem",
       x = "% Reaction Change",
       y = "subsystem") +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 10))
#ggsave(filename = "pct_changed_subsystems.pdf", plot = last_plot(), width = 8, height = 8, units = "in", dpi = 300)


