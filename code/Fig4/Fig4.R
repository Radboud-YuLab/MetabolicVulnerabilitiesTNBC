# plot heatmaps for rxns fluxes that are significantly different
# Programmer: Cyriel Huijer
# Updated: 22 June 2026

#packages
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)

# clear environment
rm(list=ls())

# set wd to /code
setwd(dir = "")

translation_table <- read.csv("data/modeling/humanGEM2BiGG.csv")

convert_ensembl_to_symbol <- function(x) {
  
  if (is.na(x)) return(NA)
  # detect original separator
  sep <- if (grepl(" and ", x)) "&" else "/"
  # split by both " or " and " and "
  ids <- unlist(strsplit(x, " and | or "))
  # map to gene symbols
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  # combine back with the  separator
  paste(symbols, collapse = sep)
}

plot_flux_heatmap <- function(df, title_name, max_value) {
  
  # Order dataframe by subsystem
  df <- df[order(df$subsystem), ]
  
  # Matrix for heatmap
  mat <- as.matrix(df[, 5:6])
  mat <- abs(mat)
  rownames(mat) <- df$bigg_rxn
  colnames(mat) <- c("0", "1000")
  
  # Rxn count as factor
  rxn_count_fac <- factor(df$rxn_count, levels = 1:4)
  
  # Row annotation
  row_ha <- rowAnnotation(
    Subsystem = df$subsystem,
    RxnCount = rxn_count_fac,
    col = list(
      Subsystem = subsystem_colors,
      RxnCount = c("1" = "#f0f0f0",
                   "2" = "#bdbdbd",
                   "3" = "#636363",
                   "4" = "#000000")
    ),
    show_annotation_name = TRUE
  )
  
  # Column color
  col_fun <- colorRamp2(
    c(0, max_value),
    c("grey90", "red")
  )
  
  # Heatmap
  Heatmap(mat,
          name = "flux",
          col = col_fun,
          cluster_rows = FALSE,   # already ordered by subsystem
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          right_annotation = row_ha,
          column_title = title_name)
}

# from extract_rs_results.mat
rxns_myc <- read.csv("data/modeling/output/alt_rxns/20260211_MYC_rxns.csv")
rxns_175 <- read.csv("data/modeling/output/alt_rxns/20260211_175_rxns.csv")
rxns_248 <- read.csv("data/modeling/output/alt_rxns/20260211_248_rxns.csv")
rxns_273 <- read.csv("data/modeling/output/alt_rxns/20260211_273_rxns.csv")

# make list, this saves redundant code
rxn_list <- list(
  MYC = rxns_myc,
  X175 = rxns_175,
  X248 = rxns_248,
  X273 = rxns_273
)

# FIlter transport reactions:
rxn_list <- lapply(rxn_list, function(df) {
  df %>% filter(!subsystem %in% c("Transport reactions","Exchange/demand reactions"))
})

#sort high to low
rxn_list <- lapply(rxn_list, function(df) {
  df %>%
    arrange(desc(abs(raw_diff)))
})

# select top 10 rxns:
rxn_list <- lapply(rxn_list, function(df) {
  df %>%
    arrange(desc(abs(raw_diff))) %>%
    slice_head(n = 10)
})

# add missing rxn names
rxn_list <- lapply(rxn_list, function(df) {
  df <- df %>%
    mutate(rxnName = case_when(
      rxn == "MAR06627" ~ "ATP:pyruvate 2-O-phosphotransferase",
      rxn == "MAR06614" ~ "ATP:nucleoside-diphosphate phosphotransferase",
      TRUE ~ rxnName  # keep existing name for all other reactions
    ))
  df
})
# clear environment
rm(rxns_myc,rxns_175,rxns_248,rxns_273)

rxns_myc_filter <- rxn_list[["MYC"]]
rxns_175_filter <- rxn_list[["X175"]]
rxns_248_filter <- rxn_list[["X248"]]
rxns_273_filter <- rxn_list[["X273"]]

# Get gene names:
rxns_myc_filter <- rxns_myc_filter %>% mutate(associated_genes_translated = sapply(associated_genes, convert_ensembl_to_symbol))
rxns_175_filter <- rxns_175_filter %>% mutate(associated_genes_translated = sapply(associated_genes, convert_ensembl_to_symbol))
rxns_248_filter <- rxns_248_filter %>% mutate(associated_genes_translated = sapply(associated_genes, convert_ensembl_to_symbol))
rxns_273_filter <- rxns_273_filter %>% mutate(associated_genes_translated = sapply(associated_genes, convert_ensembl_to_symbol))

rxns_myc_filter <- left_join(rxns_myc_filter, translation_table, by = 'rxn')
rxns_myc_filter <- rxns_myc_filter %>% mutate(bigg_rxn = ifelse(bigg_rxn == "", associated_genes_translated, bigg_rxn))
rxns_175_filter <- left_join(rxns_175_filter, translation_table, by = 'rxn')
rxns_175_filter <- rxns_175_filter %>% mutate(bigg_rxn = ifelse(bigg_rxn == "", associated_genes_translated, bigg_rxn))
rxns_248_filter <- left_join(rxns_248_filter, translation_table, by = 'rxn')
rxns_248_filter <- rxns_248_filter %>% mutate(bigg_rxn = ifelse(bigg_rxn == "", associated_genes_translated, bigg_rxn))
rxns_273_filter <- left_join(rxns_273_filter, translation_table, by = 'rxn')
rxns_273_filter <- rxns_273_filter %>% mutate(bigg_rxn = ifelse(bigg_rxn == "", associated_genes_translated, bigg_rxn))


#common_rxns <- Reduce(intersect, list(rxns_myc_filter$rxn,rxns_175_filter$rxn,rxns_248_filter$rxn,rxns_273_filter$rxn))

# Add common rxns column:
#rxns_myc_filter$common_rxn  <- rxns_myc_filter$rxn  %in% common_rxns
#rxns_175_filter$common_rxn  <- rxns_175_filter$rxn  %in% common_rxns
#rxns_248_filter$common_rxn  <- rxns_248_filter$rxn  %in% common_rxns
#rxns_273_filter$common_rxn  <- rxns_273_filter$rxn  %in% common_rxns

rxn_vectors <- list(rxns_myc_filter$rxn,rxns_175_filter$rxn,rxns_248_filter$rxn, rxns_273_filter$rxn)

# Create a named vector with counts for all unique reactions
all_rxns <- unique(unlist(rxn_vectors))

# Count in how many dataframes each rxn appears
rxn_counts <- sapply(all_rxns, function(x) sum(sapply(rxn_vectors, function(v) x %in% v)))

# Add count column to each dataframe
rxns_myc_filter$rxn_count <- rxn_counts[rxns_myc_filter$rxn]
rxns_175_filter$rxn_count <- rxn_counts[rxns_175_filter$rxn]
rxns_248_filter$rxn_count <- rxn_counts[rxns_248_filter$rxn]
rxns_273_filter$rxn_count <- rxn_counts[rxns_273_filter$rxn]

## For enzyme lookup:
rxn_list_filtered <- list(
  MYC  = rxns_myc_filter,
  X175 = rxns_175_filter,
  X248 = rxns_248_filter,
  X273 = rxns_273_filter
)

# Combine all, keep only needed columns, and get unique reactions
unique_rxns_df <- bind_rows(rxn_list_filtered) %>%
  dplyr::select(rxn, rxnName, associated_genes_translated) %>%
  distinct() # keep only unique rows
#write.csv(unique_rxns_df,"modeling/alt_rxns/rxns_to_check.csv",row.names = F)

## make heatmaps --------------------------------------------------------------

all_subsystems <- unique(c(
  rxns_myc_filter$subsystem,
  rxns_175_filter$subsystem,
  rxns_248_filter$subsystem,
  rxns_273_filter$subsystem
))
fixed_colors <- c("#D81B60","#1E88E5","#FFC107","#004D40","#673DEC","#523175","#686C96","#F366FB","#D19F62","#E0DDC5","#600B3A","#2D3EAB","#A8E3E2","#D55E0F","#CD9BEA","#28DB7B","#DCFCAA")
fixed_colors <- fixed_colors[1:length(all_subsystems)]
subsystem_colors <- setNames(fixed_colors, all_subsystems)

ht_myc  <- plot_flux_heatmap(rxns_myc_filter,"MYC",15)
ht_175  <- plot_flux_heatmap(rxns_175_filter,"X175",15)
ht_248  <- plot_flux_heatmap(rxns_248_filter,"X248",15)
ht_273  <- plot_flux_heatmap(rxns_273_filter,"X273",15)

ht_myc
ht_175
ht_248
ht_273

# Fig 4A
#pdf("figures/final_figures/Fig7/ht_myc_high.pdf",  width = 7, height = 5)
draw(ht_myc)
#dev.off()

# Fig 4B
#pdf("figures/final_figures/Fig7/ht_175_high.pdf",  width = 7, height = 5)
draw(ht_175)
#dev.off()

# Fig 4C
#pdf("figures/final_figures/Fig7/ht_248_high.pdf",  width = 7, height = 5)
draw(ht_248)
#dev.off()

# Fig 4D
#pdf("figures/final_figures/Fig7/ht_273_high.pdf",  width = 7, height = 5)
draw(ht_273)
#dev.off()

## now do it other, take log2fc (Fig 4E-F: ------------------------------------------------

rxns_myc <- read.csv("data/modeling/output/alt_rxns/20260211_MYC_rxns.csv")
rxns_175 <- read.csv("data/modeling/output/alt_rxns/20260211_175_rxns.csv")
rxns_248 <- read.csv("data/modeling/output/alt_rxns/20260211_248_rxns.csv")
rxns_273 <- read.csv("data/modeling/output/alt_rxns/20260211_273_rxns.csv")

# make list, this saves redundant code
rxn_list <- list(
  MYC = rxns_myc,
  X175 = rxns_175,
  X248 = rxns_248,
  X273 = rxns_273
)

# FIlter transport reactions:
rxn_list <- lapply(rxn_list, function(df) {
  df %>% filter(!subsystem %in% c("Transport reactions","Exchange/demand reactions"))
})

rxn_list <- lapply(rxn_list, function(df) {
  col_0 <- grep("^mean_sample_.*_0$", colnames(df), value = TRUE)
  col_1000 <- grep("^mean_sample_.*_1000$", colnames(df), value = TRUE)
  df %>% mutate(log2FC_abs = log2((abs(.data[[col_1000]])) /(abs(.data[[col_0]]))))
})

# select top 10 rxns:
rxn_list <- lapply(rxn_list, function(df) {
  df %>%
    arrange(desc(abs(log2FC_abs))) %>%
    slice_head(n = 10)
})

# add missing rxn names
rxn_list <- lapply(rxn_list, function(df) {
  df <- df %>%
    mutate(rxnName = case_when(
      rxn == "MAR06627" ~ "ATP:pyruvate 2-O-phosphotransferase",
      rxn == "MAR06614" ~ "ATP:nucleoside-diphosphate phosphotransferase",
      rxn == "MAR03819" ~ "L-glutamate gamma-semialdehyde dehydrogenase",
      rxn == "MAR08452" ~ "ATP:(d)CMP/UMP phosphotransferase",
      rxn == "MAR02633" ~ "palmitoyl-CoA:L-carnitine O-palmitoyltransferase",
      rxn == "MAR03969" ~ "ATP:uridine/cytidine 5'-phosphotransferase",
      rxn == "MAR03819" ~ "L-glutamate gamma-semialdehyde dehydrogenase",
      rxn == "MAR08450" ~ "ATP:(d)CMP/UMP phosphotransferase (MAR08450)",
      rxn == "MAR08456" ~ "ATP:(d)CMP/UMP phosphotransferase (MAR08456)",
      TRUE ~ rxnName  # keep existing name for all other reactions
    ))
  df
})
# clear environment
rm(rxns_myc,rxns_175,rxns_248,rxns_273)

rxns_myc_filter <- rxn_list[["MYC"]]
rxns_175_filter <- rxn_list[["X175"]]
rxns_248_filter <- rxn_list[["X248"]]
rxns_273_filter <- rxn_list[["X273"]]

# Get gene names:
rxns_myc_filter <- rxns_myc_filter %>% mutate(associated_genes_translated = sapply(associated_genes, convert_ensembl_to_symbol))
rxns_175_filter <- rxns_175_filter %>% mutate(associated_genes_translated = sapply(associated_genes, convert_ensembl_to_symbol))
rxns_248_filter <- rxns_248_filter %>% mutate(associated_genes_translated = sapply(associated_genes, convert_ensembl_to_symbol))
rxns_273_filter <- rxns_273_filter %>% mutate(associated_genes_translated = sapply(associated_genes, convert_ensembl_to_symbol))

rxns_myc_filter <- left_join(rxns_myc_filter, translation_table, by = 'rxn')
rxns_myc_filter <- rxns_myc_filter %>% mutate(bigg_rxn = ifelse(bigg_rxn == "", associated_genes_translated, bigg_rxn))
rxns_175_filter <- left_join(rxns_175_filter, translation_table, by = 'rxn')
rxns_175_filter <- rxns_175_filter %>% mutate(bigg_rxn = ifelse(bigg_rxn == "", associated_genes_translated, bigg_rxn))
rxns_248_filter <- left_join(rxns_248_filter, translation_table, by = 'rxn')
rxns_248_filter <- rxns_248_filter %>% mutate(bigg_rxn = ifelse(bigg_rxn == "", associated_genes_translated, bigg_rxn))
rxns_273_filter <- left_join(rxns_273_filter, translation_table, by = 'rxn')
rxns_273_filter <- rxns_273_filter %>% mutate(bigg_rxn = ifelse(bigg_rxn == "", associated_genes_translated, bigg_rxn))

rxn_vectors <- list(rxns_myc_filter$rxn,rxns_175_filter$rxn,rxns_248_filter$rxn, rxns_273_filter$rxn)

# Create a named vector with counts for all unique reactions
all_rxns <- unique(unlist(rxn_vectors))

# Count in how many dataframes each rxn appears
rxn_counts <- sapply(all_rxns, function(x) sum(sapply(rxn_vectors, function(v) x %in% v)))

# Add count column to each dataframe
rxns_myc_filter$rxn_count <- rxn_counts[rxns_myc_filter$rxn]
rxns_175_filter$rxn_count <- rxn_counts[rxns_175_filter$rxn]
rxns_248_filter$rxn_count <- rxn_counts[rxns_248_filter$rxn]
rxns_273_filter$rxn_count <- rxn_counts[rxns_273_filter$rxn]

all_subsystems <- unique(c(
  rxns_myc_filter$subsystem,
  rxns_175_filter$subsystem,
  rxns_248_filter$subsystem,
  rxns_273_filter$subsystem
))

all_subsystems <- unique(c(
  rxns_myc_filter$subsystem,
  rxns_175_filter$subsystem,
  rxns_248_filter$subsystem,
  rxns_273_filter$subsystem
))

#fixed_colors <- c("#D81B60","#1E88E5","#FFC107","#004D40","#673DEC","#523175","#686C96","#F366FB","#D19F62","#E0DDC5","#600B3A","#2D3EAB","#A8E3E2","#D55E0F","#CD9BEA","#28DB7B","#DCFCAA")
fixed_colors <- c("#004D40","#673DEC","#523175","#686C96","#F366FB","#FFC107","#1E88E5","#D19F62","#E0DDC5","#600B3A","#2D3EAB","#A8E3E2","#D55E0F","#CD9BEA","#28DB7B","#DCFCAA","#D81B60")

fixed_colors <- fixed_colors[1:length(all_subsystems)]
subsystem_colors <- setNames(fixed_colors, all_subsystems)


ht_myc  <- plot_flux_heatmap(rxns_myc_filter,"MYC",max(abs(rxns_myc_filter[, 5:6])))
ht_175  <- plot_flux_heatmap(rxns_175_filter,"X175",max(abs(rxns_175_filter[, 5:6])))
ht_248  <- plot_flux_heatmap(rxns_248_filter,"X248",max(abs(rxns_248_filter[, 5:6])))
ht_273  <- plot_flux_heatmap(rxns_273_filter,"X273",max(abs(rxns_273_filter[, 5:6])))

ht_myc
ht_175
ht_248
ht_273

# Fig 4E
#pdf("figures/final_figures/Fig7/ht_myc_low.pdf", width = 7, height = 5)
draw(ht_myc)
#dev.off()

# Fig 4F
#pdf("figures/final_figures/Fig7/ht_175_low.pdf", width = 7, height = 5)
draw(ht_175)
#dev.off()

# Fig 4G
#pdf("figures/final_figures/Fig7/ht_248_low.pdf", width = 7, height = 5)
draw(ht_248)
#dev.off()

# Fig 4H
#pdf("figures/final_figures/Fig7/ht_273_low.pdf", width = 7, height = 5)
draw(ht_273)
#dev.off()




