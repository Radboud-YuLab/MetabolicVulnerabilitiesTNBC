## exhcange fluxes in 2D plot:
# written by: Cyriel Huijer 
# Last updated: 15 June 2026

library(tidyverse)
library(purrr)
library(broom)

# Clear environment
rm(list=ls())

# Set working directory to /code
setwd("")

data <- read.csv("data/exchange_fluxes/exch_fluxes_2D.csv") # produced in 20250310_flux_determination_CH00027.R
data$cell_line[data$cell_line == "E"] <- "EV"
data$cell_line[data$cell_line == "M"] <- "MYC"
data$cell_line[data$cell_line == "T"] <- "TP53R175H"
data$cell_line[data$cell_line == "TM"] <- "TP53R175H_MYC"


data$condition_linked <- paste0(data$cell_line, "_", data$dox_conc)

## Stats ANOVA

# set factors for order of bars
data$cell_line <- factor(data$cell_line)
data$amino_acid <- factor(data$amino_acid)
data$condition_linked <- factor(data$condition_linked)

tableS1 <- data %>%
  dplyr::select(condition_linked,amino_acid,replicate, slope, R2)
colnames(tableS1) <- c("condition", "metabolite","replicate","flux_mmol_cell_h","R2")
#write.csv(tableS1, "../supplemental_tables/20260615_TableS1.csv", row.names = FALSE)

# Check anova per condition
models <- data %>% 
  group_by(amino_acid, cell_line) %>%
  group_modify(~ {
    model <- aov(slope ~ condition_linked, data = .x)
    
    tibble(model = list(model))
  }) %>%
  ungroup()

anova_result <- models %>%
  mutate(anova = purrr::map(model, broom::tidy)) %>%
  unnest(anova) %>%
  filter(term == "condition_linked")

tukey_result <- models %>%
  mutate(tukey = purrr::map(model, ~ {
    TukeyHSD(.x, "condition_linked")$condition_linked %>%
      as.data.frame() %>%
      tibble::rownames_to_column("comparison") %>%
      tidyr::separate(comparison, into = c("group1", "group2"), sep = "-")
  })) %>%
  unnest(tukey)


# THe order of bars
order_vec <- c(
  "EV_0", "EV_100","EV_1000",
  "MYC_0", "MYC_100","MYC_1000",
  "TP53R175H_0","TP53R175H_100","TP53R175H_1000",
  "TP53R175H_MYC_0","TP53R175H_MYC_100","TP53R175H_MYC_1000"
)

colors <- c(
  "MYC" = "#004D40",
  "EV"  = "#4EC6A9",  
  "TP53R175H" = "#6A0E97",
  "TP53R175H_MYC" = "#D81B60"  
)



## Plotting

# Filter for amino acid Glc
df_glc <- data %>%
  filter(amino_acid == "Glc")

# calculate means
means_slopes <- df_glc %>%
  group_by(subline = cell_line, sub_dox = condition_linked) %>%
  summarise(mean_slope = mean(slope, na.rm = TRUE), .groups = "drop")

# use factors to get correct order of bars in plot
means_slopes$sub_dox <- factor(means_slopes$sub_dox, levels = order_vec)
df_glc$condition_linked <- factor(df_glc$condition_linked, levels = order_vec)

# Plot
p_Glc <- ggplot() +
  
  # bar plot (means)
  geom_col(data = means_slopes,
           aes(x = sub_dox, y = mean_slope, fill = subline),
           width = 0.8) +
  
  # Points (individual slopes)
  geom_point(data = df_glc,
             aes(x = condition_linked, y = slope),
             size = 2,
             fill = "black",
             color = "black",
             stroke = 0.8) +
  
  scale_fill_manual(values = colors) +
  
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  
  labs(y = "Glucose consumption (mmol/cell/hr)", x = "")

all_plots <- list()

for (aa in unique(data$amino_acid)) {
  # subset data
  df_sub <- data %>% filter(amino_acid == aa)
  
  # calculate mean slopes per cell line & condition
  means_slopes <- df_sub %>%
    group_by(subline = cell_line, sub_dox = condition_linked) %>%
    summarise(mean_slope = mean(slope, na.rm = TRUE), .groups = "drop")
  
  # set factors
  means_slopes$sub_dox <- factor(means_slopes$sub_dox, levels = order_vec)
  df_sub$condition_linked <- factor(df_sub$condition_linked, levels = order_vec)
  
  # plot
  p <- ggplot() +
    geom_col(data = means_slopes,
             aes(x = sub_dox, y = mean_slope, fill = subline),
             width = 0.8) +
    geom_point(data = df_sub,
               aes(x = condition_linked, y = slope),
               size = 2, fill = "black", color = "black", stroke = 0.8) +
    scale_fill_manual(values = colors) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      y = paste0(aa, " consumption (mmol/cell/hr)"),
      x = "",
      title = aa
    )
  
  # store plot in list with amino acid name
  all_plots[[aa]] <- p
}

for (i in unique(data$amino_acid)){
  print(all_plots[[i]])
}

# Fig 2B
all_plots[["Glc"]]
#ggsave(filename = "figures/final_figures/Fig5/Glc_consumption_2D.pdf", plot = last_plot(), width = 6, height = 4, units = "in", dpi = 300)
# Fig 2C
all_plots[["Lac"]]
#ggsave(filename = "figures/final_figures/Fig5/Thr_consumption_2D.pdf", plot = last_plot(), width = 6, height = 4, units = "in", dpi = 300)








