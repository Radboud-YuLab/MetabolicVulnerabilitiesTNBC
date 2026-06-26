## Plot PCA of solutions random sampling:
# written by: Cyriel Huijer
# Last updated 22 June 2026

library(tidyverse)
library(ggrepel)

# clear environment
rm(list=ls())

# Set working directory to /code
setwd()

# Load data
pca_data <- read.csv(file = "data/modeling/output/20251217_simulated_mean_fluxes_per_condition.csv")

# transpose so columns = features, rows = samples
data_t <- t(pca_data) 
data_t <- data_t[, apply(data_t, 2, var) != 0]

# perform PCA
pca_res <- prcomp(data_t, scale. = TRUE)

# % variance explained
explained <- summary(pca_res)$importance[2, 1:2] * 100

# Create a data frame for plotting
plot_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  sample = colnames(pca_data)  # labels = original column names
)

colors <- c(
  "TP53R175H_MYC" = "#D81B60", 
  "TP53R248Q_MYC" = "#1E88E5",
  "TP53R273H_MYC" = "#FFC107",
  "MYC" = "#004D40", 
  "EV" = "#4EC6A9"
)

# create grouping variable based on sample names
plot_df$group <- NA
plot_df$group[grepl("TP53R175H", plot_df$sample)] <- "TP53R175H_MYC"
plot_df$group[grepl("TP53R248Q", plot_df$sample)] <- "TP53R248Q_MYC"
plot_df$group[grepl("TP53R273H", plot_df$sample)] <- "TP53R273H_MYC"
plot_df$group[grepl("^MYC", plot_df$sample)] <- "MYC"
plot_df$group[grepl("^EV", plot_df$sample)] <- "EV"



# Fig 3A
ggplot(plot_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = sample), color = "black", size = 4) +
  scale_color_manual(values = colors) +
  xlab(paste0("PC1 (", round(explained[1],0), "%)")) +
  ylab(paste0("PC2 (", round(explained[2],0), "%)")) +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank())
#ggsave(filename = "../figures/final_figures/Fig6/PCA_rxns.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)



