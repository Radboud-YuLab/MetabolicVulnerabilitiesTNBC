## qPCR plots TP53 and MYC
# written by: Cyriel Huijer
# lastupdated = 01/12/2025

library(tidyverse)

# Clear environemnt
rm(list=ls())

# Set to /code working directory
setwd(dir = "")

# load data:

data <- read.csv(file = "data/qPCR/qPCR_feed_R.csv")

# Create condition column
data$condition_linked <- paste0(data$Condition, "_", data$Dox)

# THe order of bars
order_vec <- c(
  "EVEV_0", "EVEV_1000",
  "MYC_0",  "MYC_1000",
  "TP53R175H_0","TP53R175H_1000", 
  "TP53R248Q_0","TP53R248Q_1000",
  paste0("TP53R175H_MYC_", c(0,80,160,1000)),
  paste0("TP53R248Q_MYC_", c(0,80,160,1000)),
  paste0("TP53R273H_MYC_", c(0,80,160,1000))
)

# Calculate means for both TP53 and MYC
means <- data %>%
  group_by(condition_linked, Condition) %>%
  summarise(
    mean_TP53 = mean(FC_TP53, na.rm = TRUE),
    mean_MYC  = mean(FC_MYC,  na.rm = TRUE),
    .groups = "drop"
  )

# Apply factor ordering to both data and means!
data$condition_linked <- factor(data$condition_linked, levels = order_vec)
means$condition_linked <- factor(means$condition_linked, levels = order_vec)

# Define your custom palette, mapping colors to unique Condition levels
colors <- c(
  "TP53R175H_MYC" = "#D81B60", 
  "TP53R248Q_MYC" = "#1E88E5", 
  "TP53R273H_MYC" = "#FFC107",
  "TP53R248Q" = "#7E4E81",
  "TP53R175H" = '#6A0E97',
  "MYC" = "#004D40", 
  "EVEV" = "#4EC6A9"
)

# Plot
plot_tp53 <- ggplot(data, aes(x = condition_linked, y = FC_TP53)) +
  
  # barplot for mean
  geom_col(data = means,
           aes(x = condition_linked, y = mean_TP53, fill = Condition),
           alpha = 1,
           width = 0.8) +
  
  # replicate points
  geom_point(size = 3, alpha = 0.9) +
  scale_fill_manual(values = colors) +
  
  labs(title = "TP53 overexpression",
       y = "Fold-change diff. to control") +
  
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_tp53

# Plot
plot_myc <- ggplot(data, aes(x = condition_linked, y = FC_MYC)) +
  
  # barplot for mean
  geom_col(data = means,
           aes(x = condition_linked, y = mean_MYC, fill = Condition),
           alpha = 1,
           width = 0.8) +
  
  # Replicate points
  geom_point(size = 3) +
  scale_fill_manual(values = colors) +
  
  labs(title = "MYC overexpression",
       y = "Fold-change diff. to control") +
  
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_myc

# save plots
# Fig 1B
#ggsave("TP53_overexpression_qpcr.pdf", plot = plot_tp53, width = 12, height = 6, units = "in", dpi = 300)

# Fig 1A
#ggsave("MYC_overexpression_qpcr.pdf", plot = plot_myc, width = 12, height = 6, units = "in", dpi = 300)
