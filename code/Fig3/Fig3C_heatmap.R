## Plot heatmap of rxns:
# Programmer: Cyriel Huijer
# Updated: 22 June 2026
library(tidyverse)
library(ComplexHeatmap)

#clean environment
rm(list=ls())

# set wd to /code
setwd("")

# load in data 
# min-max normalized
data <- read.csv(file = "data/modeling/output/20251218_mean_rxns_heatmap.csv")
data[,1:18] <- abs(data[,1:18])
sub_check <- as.data.frame(table(data$subsystem))
sub_check <- sub_check %>% filter(Freq > 5 & Freq < 21)
data <- data %>% filter(subsystem %in% sub_check$Var1)

# min-max normalize the ata
normalized_data <- data
for (i in 1:nrow(data)){
  values <- data[i,1:18]
  min_val <- min(values)
  max_val <- max(values)
  for (j in 1:length(values)){
    norm_value <- (values[j]-min_val)/(max_val - min_val)
    normalized_data[i,j] <- norm_value
  }
}

data <- normalized_data
data <- data[order(data$subsystem), ]

# extract matrix for heatmap (numeric matrix, rows x columns)
mat <- as.matrix(data[, 1:18])

# Row annotation: factor of subsystems, keeping the order
subsystem_factor <- factor(data$subsystem, levels = unique(data$subsystem))

# Generate colors for subsystems (one distinct color per level)
subsystem_colors <- structure(
  circlize::rand_color(length(levels(subsystem_factor))),
  names = levels(subsystem_factor)
)

row_ha <- rowAnnotation(
  Subsystem = subsystem_factor,
  col = list(Subsystem = subsystem_colors),
  show_annotation_name = TRUE
)

col_fun <- colorRamp2(c(0, 1), c("grey90", "red"))

# Fig 3C
heatmap <- Heatmap(
  mat,
  name = "min-max value",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = T,
  show_row_names = F,
  show_column_names = TRUE,
  top_annotation = NULL,
  right_annotation = row_ha,
  row_split = subsystem_factor
)
heatmap

#pdf("heatmap_min_max_normalized_intracellular_fluxes.pdf", width = 10, height = 5)  # dimensions in inches

draw(heatmap)

#dev.off()
