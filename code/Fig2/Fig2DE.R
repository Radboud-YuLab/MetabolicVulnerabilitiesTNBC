## 3D exchange fluxes, fig 2DE
# written by: Cyriel Huijer 
# Last updated: 22 June 2026

# packages
library(tidyverse)
library(openxlsx)
library(patchwork)

# Clear environment
rm(list=ls())

# Set working directory to /code
setwd("")

# Load data

exch_fluxes <- read.csv("data/exchange_fluxes/exch_fluxes_3D/exch_fluxes_3D.csv")

colors <- c(
  "TP53R175H_MYC" = "#D81B60", 
  "TP53R248Q_MYC" = "#1E88E5",
  "TP53R273H_MYC" = "#FFC107",
  "MYC" = "#004D40", 
  "EV" = "#4EC6A9"
)

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
plots <- exch_fluxes %>%
  filter(!condition_linked %in% skip_conditions) %>%
  split(.$metabolite) %>%
  map(~ {
    df <- .x
    
    # Use factors for correct order
    df$condition_linked <- factor(df$condition_linked, levels = order_vec)
    means_sub <- means_slope %>%
      filter(metabolite == unique(df$metabolite),
             !condition_linked %in% skip_conditions)
    means_sub$condition_linked <- factor(means_sub$condition_linked, levels = order_vec)
    
    ggplot(df, aes(x = condition_linked, y = slope)) +
      geom_col(data = means_sub, aes(x = condition_linked, y = mean_slope, fill = cell_line), width = 0.8) +
      geom_point(color = "black", fill = "black", size = 3) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      labs(
        title = unique(df$metabolite),
        y = "mmol/cell/hr",
        x = ""
      ) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })

# Fig 2D
print(plots[["Glucose"]]) # yes
# Fig 2E
print(plots[["Lactate"]]) # yes





