## 3D exchange fluxes, fig 2A
# written by: Cyriel Huijer 
# Last updated: 22 June 2026

# packages
library(tidyverse)
library(openxlsx)
library(patchwork)

# Clear environment
rm(list=ls())

# Set working directory to /code
setwd("C:/Users/chuijer/Documents/nextcloud_backup/PhD/Projects/Project1/1_manuscript/code")

# Load data

exch_fluxes <- read.csv("data/exchange_fluxes/exch_fluxes_3D/exch_fluxes_3D.csv")

# THe order of bars
order_vec <- c(
  "EV_0", "EV_1000",
  "MYC_0","MYC_80","MYC_160",  "MYC_1000",
  paste0("TP53R175H_MYC_", c(0,80,160,1000)),
  paste0("TP53R248Q_MYC_", c(0,80,160,1000)),
  paste0("TP53R273H_MYC_", c(0,80,160,1000))
)

means_slope <- exch_fluxes %>%
  group_by(condition_linked, cell_line,metabolite) %>%
  summarise(mean_slope = mean(slope, na.rm = TRUE)) %>%
  ungroup()
# this was coded in for skipping conditions, but all conditions were included eventually so "empty" now
skip_conditions <- 'empty'
# Mets that were included in the analysis after QC
yes_mets <- c("Glucose", "Lactate", "Asn",
              "Arg", "Thr", "Met",
              "Ser", "Trp", "Val")

exch_fluxes_heatmap <- means_slope %>%
  filter(!condition_linked %in% skip_conditions) %>%
  filter(metabolite %in% yes_mets)

exch_fluxes_heatmap_wide <- exch_fluxes_heatmap %>%
  dplyr::select(condition_linked, metabolite, mean_slope) %>%
  pivot_wider(names_from = metabolite, values_from = mean_slope) %>%
  as.data.frame()

row.names(exch_fluxes_heatmap_wide) <-exch_fluxes_heatmap_wide$condition_linked
exch_fluxes_heatmap_wide <- exch_fluxes_heatmap_wide %>% dplyr::select(-condition_linked)
# make matrix
mat <- as.matrix(exch_fluxes_heatmap_wide)
mat <- abs(mat)
mat_z <- scale(mat, center = TRUE, scale = TRUE)
mat_z_df <- as.data.frame(mat_z)
# mat_z_df rows are ordered according to order_vec
mat_z_df_ordered <- mat_z_df[order_vec, ]

# Convert to matrix 
mat_z <- as.matrix(mat_z_df_ordered)
mat_z_rotated <- t(mat_z_df_ordered)

col_cell_line <- c(
  rep("EV", 2),
  rep("MYC", 4),
  rep("TP53R175H_MYC", 4),
  rep("TP53R248Q_MYC", 4),
  rep("TP53R273H_MYC", 4)
)

# Name the vector with the matrix column names
names(col_cell_line) <- colnames(mat_z)

# colors for ehatmap
cell_line_colors <- c(
  "TP53R175H_MYC" = "#D81B60",
  "TP53R248Q_MYC" = "#1E88E5",
  "TP53R273H_MYC" = "#FFC107",
  "MYC" = "#004D40",
  "EV" = "#4EC6A9"
)

# annotation of columns
col_ha <- HeatmapAnnotation(
  CellLine = col_cell_line,
  col = list(CellLine = cell_line_colors),
  show_annotation_name = TRUE
)

# Plot with annotation
heatmap <- Heatmap(mat_z_rotated,
                   name = "Z-score",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   top_annotation = col_ha,
                   show_row_names = TRUE,
                   show_column_names = TRUE)
# Open PDF device with width=8, height=5 inches
#pdf("./1_manuscript/figures/final_figures/Fig5/20251216_heatmap_3D.pdf", width = 8, height = 5)

# Draw Fig 2A
draw(heatmap)




