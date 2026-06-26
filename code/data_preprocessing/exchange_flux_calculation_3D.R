## 3D exchange fluxes
# written by: Cyriel Huijer 
# Last updated: 22 June 2026

# packages
library(tidyverse)
library(openxlsx)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

# Clear environment
rm(list=ls())

# Set working directory to /code
setwd()
data <- read.csv(file = "data/growth_rate/celltiterglo_3D/output/growth_rates_filtered_reps_no_filter.csv")


# Get mean cell count at t24:
data$start_cells <- mean(data$n24)
data <- data %>%
  group_by(condition_linked) %>%
  mutate(mean_k = mean(k, na.rm = TRUE)) %>%
  ungroup()

data_mean <- data %>%
  group_by(condition_linked) %>%
  summarise(mean_k = mean(k, na.rm = TRUE),
            sd_k = sd(k,na.rm=T),.groups = "drop")
data_mean$lb <- data_mean$mean_k - data_mean$sd_k
data_mean$ub <- data_mean$mean_k + data_mean$sd_k
data_mean <- data_mean %>%mutate(lb = pmax(lb, 0))

ordered_conditions <- c(
  "EV_0", "EV_1000",
  "MYC_0", "MYC_80", "MYC_160", "MYC_1000",
  "TP53R175H_MYC_0", "TP53R175H_MYC_80", "TP53R175H_MYC_160", "TP53R175H_MYC_1000",
  "TP53R248Q_MYC_0", "TP53R248Q_MYC_80", "TP53R248Q_MYC_160", "TP53R248Q_MYC_1000",
  "TP53R273H_MYC_0", "TP53R273H_MYC_80", "TP53R273H_MYC_160", "TP53R273H_MYC_1000"
)
data_mean <- data_mean %>%mutate(condition_linked = factor(condition_linked, levels = ordered_conditions)) %>%arrange(condition_linked)

# Growth rates used for modeling:
#write.csv(data_mean, "data/growth_rate/celltiterglo_3/20251217_growth_rate.csv")

data <- data %>%
  group_by(condition_linked) %>%
  summarise(
    n24       = mean(n24, na.rm = TRUE),
    start_cells = mean(start_cells, na.rm = TRUE),
    mean_k    = mean(mean_k, na.rm = TRUE),
    .groups = "drop"
  )

# timepoints in metabolomics data:
tp <- c(24, 30, 36, 48, 54, 60, 72)

predicted_cells <- data %>%
  rowwise() %>%
  mutate(
    pred = list(
      tibble(
        time_point = tp,
        logN = n24 + mean_k * (tp - 24),  # ln(cell count) at time tp
        cells = exp(logN) # transform to non transformed cell count
      )
    )
  ) %>%
  unnest(pred) %>%
  ungroup()

# predicted_cells will be used for the flux calculation
                                                                         
glucose<- read.xlsx("data/metabolomics/HPLC/20250603_HPLC_R1.xlsx")
lactate<- read.xlsx("data/metabolomics/HPLC/20250603_HPLC_R1.xlsx",sheet = "lactate")
hplc_correction <- glucose %>%
  filter(cell_line == "B1") %>%
  summarise(corr_factor = 3.151/mean(concentration)) %>% # correct for 
  pull(corr_factor)

glucose$condition_linked <- paste0(glucose$cell_line,"_",glucose$dox)
lactate$condition_linked <- paste0(lactate$cell_line,"_",lactate$dox)

glucose <- glucose %>%
  distinct(sample_ID, .keep_all = TRUE)
lactate <- lactate %>%
  distinct(sample_ID, .keep_all = TRUE)

glucose$rep <- 1
lactate$rep <- 1

glucose <- glucose %>% 
  filter(! cell_line %in% c("B1","MQ","STD1","STD2","STD3","STD4","STD5","STD6","STD7","STD8")) %>%
  dplyr::select(cell_line, condition_linked,rep,time_point, dox, concentration)
lactate <- lactate %>% 
  filter(! cell_line %in% c("B1","MQ","STD1","STD2","STD3","STD4","STD5","STD6","STD7","STD8")) %>%
  dplyr::select(cell_line, condition_linked,rep,time_point, dox, concentration)

# correct concentrations
glucose$corrected_concentration <- glucose$concentration * hplc_correction
lactate$corrected_concentration <- lactate$concentration * hplc_correction

# add metabolite name:
glucose$metabolite <- "Glucose"
lactate$metabolite <- "Lactate"

# mol weight glucose = 180.16 g/mol and lactate = 88.06
glucose$corrected_mmol_liter <- (glucose$corrected_concentration / 180.16) * 1000
lactate$corrected_mmol_liter <- (lactate$corrected_concentration / 88.06) * 1000

hplc_data <- rbind(glucose,lactate)
hplc_data <- hplc_data %>%
  dplyr::select(cell_line, condition_linked,rep,time_point,metabolite,corrected_mmol_liter)

## Amino acids 

amino_acids <- read.xlsx(xlsxFile = "data/metabolomics/LCMS/clean_LCMS.xlsx", sheet = "all_cell_lines")
colnames(amino_acids)[5:24] <- tools::toTitleCase(colnames(amino_acids)[5:24])

# pivot the dataframe longer to match format hplc data

amino_acids <- amino_acids %>%
  pivot_longer(
    cols = Asn:Trp,
    names_to = "metabolite",
    values_to = "corrected_umol_liter"
  )
amino_acids$corrected_mmol_liter <- amino_acids$corrected_umol_liter / 1000 
amino_acids$condition_linked <- paste0(amino_acids$cell_line,"_",amino_acids$dox)
amino_acids <- amino_acids %>%
  dplyr::select(cell_line,condition_linked,rep,time_point,metabolite,corrected_mmol_liter)

# Combine all mets
all_mets <- rbind(hplc_data, amino_acids)

## Exchange flux calculation --------------------------------------------------

all_mets$mmol <- all_mets$corrected_mmol_liter * 0.0002 # volume correction in which cells were grown to go from mmol/L to mmol (assumes 200 uL media)

# add number of cells:
all_mets <- all_mets %>%
  left_join(
    predicted_cells %>%
      dplyr::select(condition_linked, time_point, cells),
    by = c("condition_linked", "time_point")
  )

all_mets$mmol_cell <- all_mets$mmol / all_mets$cells

# Metabolites that are filtered out after manual QC: 
# Asp, is depleted after 48 hours, leading to an increased consumption of Glu
# Cells are not in steady state during this transformation, therefore, models should be run before and after
# Furhter amino acids that should be removed: Cys, Leu, Ile, Lys, Gly

all_mets <- all_mets %>%
  filter(! metabolite %in% c("Cys","Leu","Ile", "Lys","Gly"))
all_mets$filter_id <- paste0(all_mets$condition_linked,"_",all_mets$metabolite,"_",all_mets$rep,"_",all_mets$time_point)
all_mets <- all_mets %>% filter(!filter_id %in% c("TP53R175H_MYC_0_Val_2_36","EV_1000_Tyr_1_24","EV_1000_Tyr_2_24","TP53R273H_MYC_1000_Ser_1_48","TP53R273H_MYC_1000_Ser_2_48","EV_0_Met_2_72","EV_0_Met_2_48","EV_0_Met_1_48"))
all_mets <- all_mets %>%dplyr::mutate(concentration = as.numeric(stringr::str_extract(condition_linked, "[^_]+$")),concentration_f = factor(concentration,levels = c(0, 80, 160, 1000)))
conc_colors <- c(
  "0"    = "#666666",
  "80"   = "#1b9e77",
  "160"  = "#d95f02",
  "1000" = "#7570b3"
)


# Calculte excange fluxes here, it is the slope of the mmol/cell curve:

exch_fluxes <- all_mets %>%
  group_by(cell_line,condition_linked, metabolite, rep) %>%
  do({
    fit <- lm(mmol_cell ~ time_point, data = .)
    tibble(
      slope = coef(fit)[["time_point"]],
      r2 = summary(fit)$r.squared
    )
  }) %>%
  ungroup()

#Did not pass QC (Gln, Ala, Pro), the R2 is too much spread out
exch_fluxes <- exch_fluxes %>%
  filter(! metabolite %in% c("Gln","Ala","Pro"))

exch_fluxes <- exch_fluxes %>%
  filter(r2 > 0.5 | metabolite == "Lactate")

# Output exch. fluxes for fig2A, Fig 2D-E
#write.csv(exch_fluxes,"data/exchange_fluxes/exch_fluxes_3D/exch_fluxes_3D.csv", row.names = F)

means_slope <- exch_fluxes %>%
  group_by(condition_linked, cell_line,metabolite) %>%
  summarise(mean_slope = mean(slope, na.rm = TRUE)) %>%
  ungroup()

skip_conditions <- 'empty'

yes_mets <- c("Glucose", "Lactate", "Asn",
              "Arg", "Thr", "Met",
              "Ser", "Trp", "Val")


## calculation of exchange fluxes used for modeling: ----------------------------------------------------------------------
# Check anova per condition
models <- exch_fluxes %>% 
  group_by(metabolite, cell_line) %>%
  group_modify(~ {
    #check how many unique levels in this group
    n_levels <- n_distinct(.x$condition_linked)
    if (n_levels >= 2) {
      model <- aov(slope ~ condition_linked, data = .x)
      tibble(model = list(model))
    } else {
      # return empty list
      tibble(model = list(NA))
    }
  }) %>%
  ungroup()

anova_result <- models %>%
  mutate(anova = purrr::map(model, broom::tidy)) %>%
  unnest(anova) %>%
  filter(term == "condition_linked")
anova_result <- anova_result %>% dplyr::select(-model)

tukey_result <- models %>%
  filter(!is.na(model)) %>% # Remove NA models
  mutate(tukey = purrr::map(model, ~ {
    TukeyHSD(.x, "condition_linked")$condition_linked %>%
      as.data.frame() %>%
      tibble::rownames_to_column("comparison") %>%
      tidyr::separate(comparison, into = c("group1", "group2"), sep = "-")
  })) %>%
  unnest(tukey)
tukey_result <- tukey_result %>% dplyr::select(-model)

## calculate exchange fluxes:
# correct for dry weight of MCF10A cells
exch_fluxes <- exch_fluxes %>% mutate(mmol_gDW = slope * (1 / (397 * 10^-12)))

exch_fluxes <- exch_fluxes %>% filter(metabolite %in% yes_mets)

exch_fluxes <- exch_fluxes %>%
  left_join(anova_result,
            by = c("cell_line", "metabolite"))

exch_fluxes <- exch_fluxes %>% dplyr::select(cell_line, condition_linked,metabolite, rep, mmol_gDW, p.value)
# check whether a metabolite is significant, and whether its glucose or lactate (these have only one rep)
exch_fluxes_final <- exch_fluxes %>%
  mutate(
    is_special = metabolite %in% c("Glucose", "Lactate"),
    significant = p.value < 0.05
  )
# NAs at significant should be FALSE 
exch_fluxes_final <- exch_fluxes_final %>%
  mutate(
    significant = if_else(is.na(p.value), FALSE, p.value < 0.05)
  )
# if glucose or lactate, calculate mean and sd as 10% cutoff
exch_fluxes_final <- exch_fluxes_final %>%
  mutate(
    mean_final = if_else(is_special, mmol_gDW, NA_real_),
    sd_final   = if_else(is_special, abs(mean_final) * 0.10, NA_real_)
  )
# For significant values, calculate the mean and sd for condition_linked
exch_fluxes_final <- exch_fluxes_final %>%
  group_by(cell_line, metabolite, condition_linked) %>%
  mutate(
    mean_final = if_else(!is_special & significant, mean(mmol_gDW, na.rm = TRUE), mean_final),
    sd_final   = if_else(!is_special & significant, sd(mmol_gDW, na.rm = TRUE), sd_final)
  ) %>%
  ungroup()
# Now if its non signfiicnat, calculate the mean + sd on the cell_line level
exch_fluxes_final <- exch_fluxes_final %>%
  group_by(cell_line, metabolite) %>%
  mutate(
    mean_final = if_else(!is_special & !significant, mean(mmol_gDW, na.rm = TRUE), mean_final),
    sd_final   = if_else(!is_special & !significant, sd(mmol_gDW, na.rm = TRUE), sd_final)
  ) %>%
  ungroup()

# only one sd is left NA now, use the mean of all other sds for that cell line + metabolite
exch_fluxes_final <- exch_fluxes_final %>%
  group_by(cell_line, metabolite) %>%
  mutate(
    sd_final = if_else(
      is.na(sd_final),
      mean(sd_final, na.rm = TRUE),# mean of other SDs in the group
      sd_final
    )
  ) %>%
  ungroup()

exch_fluxes_summary <- exch_fluxes_final %>%
  group_by(cell_line, metabolite, condition_linked) %>%
  summarise(
    mean_final = dplyr::first(mean_final),
    sd_final   = dplyr::first(sd_final),
    .groups = "drop"
  )
# Calculate lb and ub (mean + sd and mean - sd)
exch_fluxes_summary <- exch_fluxes_summary %>%
  mutate(
    lb = mean_final - sd_final,
    ub = mean_final + sd_final
  )

# Lower bound dataframe
lb_fluxes <- exch_fluxes_summary %>%
  dplyr::select(metabolite, condition_linked, lb) %>%
  pivot_wider(
    names_from = condition_linked,
    values_from = lb
  )

# Upper bound dataframe
ub_fluxes <- exch_fluxes_summary %>%
  dplyr::select(metabolite, condition_linked, ub) %>%
  pivot_wider(
    names_from = condition_linked,
    values_from = ub
  )

# order of conditions
ordered_conditions <- c(
  "EV_0", "EV_1000",
  "MYC_0", "MYC_80", "MYC_160", "MYC_1000",
  "TP53R175H_MYC_0", "TP53R175H_MYC_80", "TP53R175H_MYC_160", "TP53R175H_MYC_1000",
  "TP53R248Q_MYC_0", "TP53R248Q_MYC_80", "TP53R248Q_MYC_160", "TP53R248Q_MYC_1000",
  "TP53R273H_MYC_0", "TP53R273H_MYC_80", "TP53R273H_MYC_160", "TP53R273H_MYC_1000"
)

# make sure 'metabolite' is first
ordered_columns <- c("metabolite", ordered_conditions)

# reorder lb_fluxes
lb_fluxes <- lb_fluxes[, ordered_columns]
# change to human gem identifiers
lb_fluxes$metabolite <- c('arginine', 'asparagine', 'glucose', 'L_lactate','methionine', 'serine', 'threonine', 'tryptophan', 'valine')
# reorder ub_fluxes
ub_fluxes <- ub_fluxes[, ordered_columns]
ub_fluxes$metabolite <- c('arginine', 'asparagine', 'glucose', 'L_lactate','methionine', 'serine', 'threonine', 'tryptophan', 'valine')

# Output for modeling
#write.csv(lb_fluxes, file = "data/exchange_fluxes/exch_fluxes_3D/for_modeling/20251217_lb.csv",row.names = F)
#write.csv(ub_fluxes, file = "data/exchange_fluxes/exch_fluxes_3D/for_modeling/20251217_ub.csv",row.names = F)




