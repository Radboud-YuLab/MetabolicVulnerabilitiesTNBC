## Calculation of exchange fluxes MCF10A data. 
# Written by: Cyriel Huijer
# Last updated: 10 March 2025

library(tidyverse)
library(ggplot2)
library(MASS)
library(robustbase)
# set to /code
setwd("")

# laod data

metab_concentration <- read.csv("data/metabolomics/2D/calculated_metabolite_concentrations_outliers_included_CH00027.csv",header = TRUE)
manual_inspection <- read.csv("data/metabolomics/2D/manual_inspection_curves.csv")
media_concentration <- read.csv("data/metabolomics/2D/DMEM_F12_media_composition_mM.csv")
metab_concentration <- metab_concentration %>%
  mutate(
    final_concentration_uM = case_when(
      amino_acid == "Glc" ~ (final_concentration * 1000) / 180.16,  # Glucose conversion was already done in previous script by 1000x multiplication
      amino_acid == "Lac" ~ (final_concentration * 1e6) / 90.08,   # Lactate conversion
      TRUE ~ final_concentration  # Keep original if not Glc or Lac
    )
  )

merged_df <- merge(metab_concentration, manual_inspection, by = c("batch", "amino_acid"))
filtered_df <- subset(merged_df, keep_manual_inspection != "n")
metab_concentration <- filtered_df

# Calculate the mean for each amino_acid and timepoint
mean_concentration <- metab_concentration %>%
  group_by(amino_acid, timepoint, cell_line, dox_conc) %>%
  summarise(mean_final_concentration_uM = mean(final_concentration_uM, na.rm = TRUE), .groups = "drop")

#no significant difference in cell count so one value is taken
cell_count <- read.csv("data/metabolomics/2D/cell_numbers_all_replicates.csv")

# cells were grown in 1.5 mL, calculate amount of metabolite in total 
metab_concentration$final_concentration_mM <- metab_concentration$final_concentration_uM/1000
metab_concentration$volume_corr <- metab_concentration$final_concentration_mM *0.0015
# cells were grown in 6-well plate, correct for this:
cell_count$cell_number_6_well <- cell_count$Cell_Number * (9.6/0.32)
merged_df <- left_join(metab_concentration, cell_count, by = c("timepoint"))
merged_df$mmol_cell <- merged_df$volume_corr / merged_df$cell_number_6_well

merged_df <- merged_df %>%
  filter(cell_line %in% c("E","M","T","TM"),dox_conc %in% c(0, 100,1000))
merged_df$condition <- paste0(merged_df$cell_line,"_",as.character(merged_df$dox_conc))

calculate_slope <- function(data) {
  fit <- lm(mmol_cell ~ timepoint, data = data)
  slope <- coef(fit)["timepoint"]
  R2 <- summary(fit)$r.squared
  return(list(slope = slope, R2 = R2))
}

slopes_df <- merged_df %>%
  group_by(condition, amino_acid, replicate) %>%
  summarize(
    slope = calculate_slope(cur_data())$slope,
    R2 = calculate_slope(cur_data())$R2,
    .groups = 'drop'
  ) %>%
  filter(R2 > 0.7)

slopes_df$cell_line <- gsub("_.*", "", slopes_df$condition)
slopes_df$dox_conc <- sub(".*_", "", slopes_df$condition)

#write_csv(slopes_df,"data/exchange_fluxes/exch_fluxes_2D.csv")







