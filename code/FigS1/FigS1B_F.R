# Supplemental figure RNAseq (Fig S1B-F)
# written by: cyriel Huijer
# last updated: 22 June 2026

library(DESeq2)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(data.table)
library(openxlsx)
library(clusterProfiler)
library(ggrepel)
library(ggbreak)
library(ggvenn)
# Clear environment
rm(list=ls())

# Set working directory to /code
setwd("")

add_genesymbol = function(genelist,decimal = TRUE){
  #Function to set ensemble names to hgnc symbol
  ##If isoform is annotated, use decimal == TRUE
  if(decimal){
    ensemble_short = as.data.frame(genelist) %>% separate(genelist,".")
  }else{
    ensemble_short = as.data.frame(genelist)
  }
  ensemble_short = AnnotationDbi::select(org.Hs.eg.db,
                                         keys=ensemble_short[,1], columns=c("SYMBOL"), keytype="ENSEMBL",multiVals = "first")
  ensemble_short = distinct(ensemble_short, ENSEMBL, .keep_all = TRUE)
  return(ensemble_short$SYMBOL)
}


#Read data + formatting
counts = read.table("data/transcriptomics/counts/experiment_2/GRCh38.p14-counts.tsv",header = TRUE)
rownames(counts) = counts$gene
colnames(counts) <- sub("^[^_]+_", "", colnames(counts))
colnames(counts) <- c("gene","TM_0_1","TM_1000_1","MYC_0_1","MYC_1000_1","TM_0_2","TM_1000_2","MYC_0_2","MYC_1000_2","TM_0_3","TM_1000_3","MYC_0_3","MYC_1000_3")
counts

counts$gene = NULL
listmeta = as.data.frame(colnames(counts))
#Generate metadata file
df = as.data.frame(strsplit(listmeta[,1],split  = "_")[[1]])
df = as.data.frame(t(df))
for(i in 1:12){
  if(i == 1){
    next
  }
  suppressWarnings({
    df = rbind(df,strsplit(listmeta[,1],split  = "_")[[i]])
  })
}
df
colnames(df) = c("gene","conc","rep")
coldata = df
coldata
rownames(coldata) = listmeta[,1]

# without second replicate:
#coldata = coldata[c(1,2,3,4,9,10,11,12),] # add 18 for also last rep MYC1000 sample 33
#counts = counts[,c(10,11,13,14,16,17)]

# MYC only:
# coldata = coldata[c(3,4,7,8,11,12),]
# counts = counts[,c(3,4,7,8,11,12)]
# dds <- DESeqDataSetFromMatrix(countData = counts,
#                               colData = coldata,
#                               design = ~conc)
# TM only:
#coldata = coldata[c(1,2,5,6,9,10),]
#counts = counts[,c(1,2,5,6,9,10)]
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~conc +gene)
# MYC vs TM overexpression:
#coldata = coldata[c(2,4,6,8,10,12),]
#counts = counts[,c(2,4,6,8,10,12)]
#dds <- DESeqDataSetFromMatrix(countData = counts,
#                              colData = coldata,
#                              design = ~gene)

dds <- DESeq(dds)
resultsNames(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("gene", "conc"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p1 <- ggplot(pcaData, aes(PC1, PC2, color=conc, shape=gene)) +
  geom_point(size=3) +
  geom_text_repel(aes(label=name)) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 15)

# Fig S1B
p1
#ggsave(filename = "../figures/final_figures/Fig2/PCA_MYC_TM.pdf", plot = p1, width = 8, height = 6)

## Analysis for MYC (Fig S1E)-----------------------------------------------------
# keep only MYC samples
myc_samples <- rownames(coldata[coldata$gene == "MYC", ])
counts_myc <- counts[, myc_samples]
coldata_myc <- coldata[myc_samples, ]

dds_myc <- DESeqDataSetFromMatrix(
  countData = counts_myc,
  colData = coldata_myc,
  design = ~ conc
)

dds_myc <- DESeq(dds_myc)

res_myc <- results(dds_myc, contrast = c("conc", "1000", "0"))
res_myc <- as.data.frame(res_myc)
res_myc$gene_symbol <- add_genesymbol(rownames(res_myc))

# Ranking score
res_myc$rank_metric <- abs(res_myc$log2FoldChange) * -log10(res_myc$padj)

# Significant UP + FC threshold
sig_up_myc <- res_myc[res_myc$padj < 0.05 & res_myc$log2FoldChange > 1, ]

# Significant DOWN + FC threshold
sig_down_myc <- res_myc[res_myc$padj < 0.05 & res_myc$log2FoldChange < -1, ]

# Top 25 ranked
top25_up_myc <- sig_up_myc[order(-sig_up_myc$rank_metric), ][1:25, ]
top25_down_myc <- sig_down_myc[order(-sig_down_myc$rank_metric), ][1:25, ]

# Correct labeling using gene_symbol
res_myc$label <- NA
res_myc$label[res_myc$gene_symbol %in% top25_up_myc$gene_symbol] <-
  res_myc$gene_symbol[res_myc$gene_symbol %in% top25_up_myc$gene_symbol]

res_myc$label[res_myc$gene_symbol %in% top25_down_myc$gene_symbol] <-
  res_myc$gene_symbol[res_myc$gene_symbol %in% top25_down_myc$gene_symbol]
# grey for NS, black for significant
res_myc$dot_color <- ifelse(
  abs(res_myc$log2FoldChange) > 1 & res_myc$padj < 0.05,
  "black",
  "lightgrey"
)
# Fig S1E
p2 <- ggplot(res_myc, aes(x = log2FoldChange, y = -log10(padj), color = dot_color)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(aes(label = label), max.overlaps = Inf, size = 5) +
  scale_color_identity() +
  theme_bw(base_size = 15) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 15))

p2
#ggsave(filename = "../figures/final_figures/Fig2/20251219_MYC_volcano.pdf", plot = p2, width = 8, height = 6)

up_myc <- res_myc %>%
  filter(padj < 0.05, log2FoldChange > 1) %>%
  pull(gene_symbol) %>%
  na.omit()

down_myc <- res_myc %>%
  filter(padj < 0.05, log2FoldChange < -1) %>%
  pull(gene_symbol) %>%
  na.omit()
## Analysis for TM (Fig S1F) ----------------------------------------------------

# keep only TM samples
tm_samples <- rownames(coldata[coldata$gene == "TM", ])
counts_tm <- counts[, tm_samples]
coldata_tm <- coldata[tm_samples, ]

dds_tm <- DESeqDataSetFromMatrix(
  countData = counts_tm,
  colData = coldata_tm,
  design = ~ conc
)

dds_tm <- DESeq(dds_tm)

res_tm <- results(dds_tm, contrast = c("conc", "1000", "0"))
res_tm <- as.data.frame(res_tm)
res_tm$gene_symbol <- add_genesymbol(rownames(res_tm))


# Ranking metric = |log2FC| * -log10(padj)
res_tm$rank_metric <- abs(res_tm$log2FoldChange) * -log10(res_tm$padj)

# Significant up/down
sig_up <- res_tm[res_tm$padj < 0.05 & res_tm$log2FoldChange > 1, ]
sig_down <- res_tm[res_tm$padj < 0.05 & res_tm$log2FoldChange < -1, ]

# Top 25 by rank metric
top25_up <- sig_up[order(-sig_up$rank_metric), ][1:25, ]
top25_down <- sig_down[order(-sig_down$rank_metric), ][1:25, ]

# Correct labeling using gene_symbol
res_tm$label <- NA
res_tm$label[res_tm$gene_symbol %in% top25_up$gene_symbol] <- res_tm$gene_symbol[res_tm$gene_symbol %in% top25_up$gene_symbol]
res_tm$label[res_tm$gene_symbol %in% top25_down$gene_symbol] <- res_tm$gene_symbol[res_tm$gene_symbol %in% top25_down$gene_symbol]

# Volcano colors: grey = NS, black = significant
res_tm$dot_color <- ifelse(
  abs(res_tm$log2FoldChange) > 1 & res_tm$padj < 0.05,
  "black",
  "lightgrey"
)

# Fig 1SF
p3 <- ggplot(res_tm, aes(x = log2FoldChange, y = -log10(padj), color = dot_color)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(aes(label = label), max.overlaps = Inf, size = 5) +
  scale_color_identity() + 
  scale_alpha(range = c(0.3, 1), guide = 'none') +
  theme_bw(base_size = 15) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 20))
p3


#ggsave(filename = "../figures/final_figures/Fig2/20251219_TM_volcano.pdf", plot = p3, width = 8, height = 6)

# Extract significantly altered genes
up_tm <- res_tm %>%
  filter(padj < 0.05, log2FoldChange > 1) %>%
  pull(gene_symbol) %>%
  na.omit()

down_tm <- res_tm %>%
  filter(padj < 0.05, log2FoldChange < -1) %>%
  pull(gene_symbol) %>%
  na.omit()

## Venn diagrams for common upregulated and downregulated genes ----------------
# Fig S1C
gene_up <- list(
  "Up MYC" = up_myc,
  "Up TM" = up_tm
)
fit <- euler(gene_up)
pdf(file = "../figures/final_figures/Fig2/up_venn.pdf", width = 6, height = 6)
plot(fit, fills = list(fill = c("#1E88E5", "#D81B60"), alpha = 0.5),
     labels = TRUE, edges = TRUE)
dev.off()
ggvenn(gene_up, fill_color = c("#1E88E5", "#D81B60"))

# fig S1D
gene_down <- list(
  "Down MYC" = down_myc,
  "Down TM" = down_tm
)
fit <- euler(gene_down)
pdf(file = "../figures/final_figures/Fig2/down_venn.pdf", width = 6, height = 6)
plot(fit, fills = list(fill = c("#1E88E5", "#D81B60"), alpha = 0.5),
     labels = TRUE, edges = TRUE)
dev.off()

ggvenn(gene_down, fill_color = c("#1E88E5", "#D81B60"))


unique_to_tm <- setdiff(up_tm, up_myc)



## GO term analysis 


# 2. Run GO enrichment
GO_corr_results <- enrichGO(
  gene         = up_tm,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# Convert results to dataframe
GO_corr_df <- as.data.frame(GO_corr_results)

# 3. Plot top GO terms
top_GO_corr <- GO_corr_df %>%
  slice_head(n = 10) %>%
  mutate(
    # Convert GeneRatio to numeric fraction
    GeneRatio_num = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    # Rank order by GeneRatio (largest at top)
    Description = factor(Description, levels = Description[order(GeneRatio_num, decreasing = FALSE)])
  )



# Bubble plot
ggplot(top_GO_corr, aes(x = Description, y = GeneRatio_num)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  scale_size_continuous(name = "Gene Count") +
  coord_flip() +
  labs(
    x = "",
    y = "Gene Ratio",
    title = "GO Enrichment: High RNA–Protein Correlating Genes (Gln Lim)",
    subtitle = "Biological Process (BP)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

























## MYC targets included: -------------------------------------------------------
# MYC targets: https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_MYC_TARGETS_V2
myc_targets <- c("AIMP2","BYSL","CBX3","CDK4","DCTPP1","DDX18","DUSP2","EXOSC5","FARSA","GNL3","GRWD1","HK2","HSPD1","HSPE1","IMP4","IPO4","LAS1L","MAP3K6","MCM4","MCM5","MPHOSPH10","MRTO4","MYBBP1A","MYC","NDUFAF4","NIP7","NOC4L","NOLC1","NOP16","NOP2","NOP56","NPM1","PA2G4","PES1","PHB1","PLK1","PLK4","PPAN","PPRC1","PRMT3","PUS1","RABEPK","RCL1","RRP12","RRP9","SLC19A1","SLC29A2","SORD","SRM","SUPV3L1","TBRG4","TCOF1","TFB2M","TMEM97","UNG","UTP20","WDR43","WDR74")
res_tm$dot_color <- ifelse(
  res_tm$gene_symbol %in% myc_targets, "red",
  ifelse(
    abs(res_tm$log2FoldChange) > 1 & res_tm$padj < 0.05,
    "black",
    "lightgrey"
  )
)

# Define alpha (opacity): significant and MYC targets opaque, others more transparent
res_tm$dot_alpha <- ifelse(
  res_tm$gene_symbol %in% myc_targets | (abs(res_tm$log2FoldChange) > 1 & res_tm$padj < 0.05),
  1,# fully opaque
  0.1 # more transparent
)

p3 <- ggplot(res_tm, aes(x = log2FoldChange, y = -log10(padj), color = dot_color)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3) +
  scale_color_identity() + 
  scale_alpha(range = c(0.3, 1), guide = 'none') +  # keep alphas as-is, no legend
  theme_bw(base_size = 15) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 20))
p3
p3 <- ggplot() +
  geom_point(data = res_tm[!(res_tm$gene_symbol %in% myc_targets), ], 
             aes(x = log2FoldChange, y = -log10(padj)), 
             color = "lightgrey", alpha = 0.2) +
  
  # Significant non-MYC genes in black
  geom_point(data = res_tm[
    res_tm$gene_symbol %in% myc_targets == FALSE & 
      abs(res_tm$log2FoldChange) > 1 & res_tm$padj < 0.05, ],
    aes(x = log2FoldChange, y = -log10(padj)),
    color = "black", alpha = 1) +
  
  # MYC targets on top, bigger and bright red
  geom_point(data = res_tm[res_tm$gene_symbol %in% myc_targets, ], 
             aes(x = log2FoldChange, y = -log10(padj)),
             color = "red", alpha = 1) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  
  # Label only top 25 genes here
  geom_text_repel(data = res_tm[!is.na(res_tm$label), ],
                  aes(x = log2FoldChange, y = -log10(padj), label = label),
                  size = 3, max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "MYC targets" = "red",
      "Significant" = "black",
      "Not significant" = "lightgrey"
    ),
    name = "Gene group"
  )+
  
  theme_bw(base_size = 15) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 20))

p3



# keep only TM samples
tm_samples <- rownames(coldata[coldata$gene == "TM", ])
counts_tm <- counts[, tm_samples]
coldata_tm <- coldata[tm_samples, ]

dds_tm <- DESeqDataSetFromMatrix(
  countData = counts_tm,
  colData = coldata_tm,
  design = ~ conc
)

dds_tm <- DESeq(dds_tm)

res_tm <- results(dds_tm, contrast = c("conc", "1000", "0"))
res_tm <- as.data.frame(res_tm)
res_tm$gene_symbol <- add_genesymbol(rownames(res_tm))

res_tm$label <- ifelse( abs(res_tm$log2FoldChange) > 1 & res_tm$padj < 0.05, res_tm$gene_symbol, NA )


# Volcano colors: grey = NS, black = significant
res_tm$dot_color <- ifelse(
  abs(res_tm$log2FoldChange) > 1 & res_tm$padj < 0.05,
  "black",
  "lightgrey"
)

p3 <-ggplot(res_tm, aes(x = log2FoldChange, y = -log10(padj), color = dot_color)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3) +
  scale_color_identity() +        # Use colors as-is
  theme_bw(base_size = 15) +
  coord_cartesian(xlim = c(-3, 3),ylim = c(0,20))
p3
ggsave(filename = "../figures/final_figures/Fig2/TM_volcano.pdf", plot = p3, width = 10, height = 6)

# Extract significantly regulated genes
up_tm <- res_tm[
  res_tm$padj < 0.05 & res_tm$log2FoldChange > 1,
  "gene_symbol"
]

down_tm <- res_tm[
  res_tm$padj < 0.05 & res_tm$log2FoldChange < -1,
  "gene_symbol"
]






#dev.off()
#Differential testing + formatting
res = results(dds, contrast=c("conc","1000","0"))
res = as.data.frame(res)
res$gene = add_genesymbol(rownames(res))
res = as.data.frame(res)
res = res[!is.na(res$padj),]
res = res[order(res$padj),]
print(paste(length(res[res$padj < 0.05 & res$log2FoldChange > 1,]$gene), "upregulated genes"))
upregulated_tm <- res[res$padj < 0.05 & res$log2FoldChange > 1,]$gene
print(paste(length(res[res$padj < 0.05 & res$log2FoldChange < -1,]$gene), "downregulated genes"))
downregulated_tm <- res[res$padj < 0.05 & res$log2FoldChange < -1,]$gene
#res[res$gene %in% c(“IL1A”,“IL1B”,“IL1RN”,“AREG”,“CD74”,“GADD45B”,“BHLHE40”,“KDM6B”,“NFKB1”,“JUN”,“B2M”,“CASP4",“NLRP3”,“HLA-DRB1"),]
options(repr.plot.width=6,repr.plot.height=8,repr.plot.res = 300)
#Visualise ageing and trained immunity genes
de_mat = assay(vsd)
rownames(de_mat)= add_genesymbol(rownames(de_mat))
de_mat = as.data.frame(de_mat)
de_mat = de_mat[rownames(de_mat) %in% c(res[res$padj < 0.05 & res$log2FoldChange > 1.2  | res$log2FoldChange < -1.2,]$gene),]
anno = coldata[,c(1,2)]
p1 = pheatmap(t(scale(t(de_mat))),show_rownames = T,show_colnames = F,annotation_col =anno,cluster_cols = TRUE)
p2 = pheatmap(t(scale(t(de_mat))),show_rownames = T,show_colnames = F,annotation_col =anno,cluster_cols = FALSE)
#pdf("myc_heatmap1.pdf",width = 6,height = 20)
p1
#dev.off()
#pdf("../tmp/myc_heatmap2.pdf”,width = 6,height = 8")
p2
#dev.off()

# Volcano Plot
res$label <- ifelse((res$log2FoldChange > 1 | res$log2FoldChange < -1) & res$padj < 0.05, res$gene, NA)
res$dot_color <- ifelse((res$log2FoldChange > 1 | res$log2FoldChange < -1) & res$padj < 0.05, "lightgrey", "black")

p_threshold <- -log10(0.05)



ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = dot_color)) +
  geom_point(alpha = 0.6) +
  geom_text_repel(aes(label = label), max.overlaps = Inf, size = 2.5) + # Add labels for significant genes
  geom_hline(yintercept = p_threshold, linetype = "dashed", color = "black") +  # Dashed line for p < 0.05
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +   # Dashed lines for log2FoldChange thresholds
  scale_color_manual(values = c("lightgrey", "black")) +  # Manual color scale
  theme_bw(base_size = 15)


# GO term analysis

res <- results(dds, contrast=c("conc","1000","0"))
res <- as.data.frame(res)
res$gene <- rownames(res)
res$gene_name <- add_genesymbol(rownames(res))

res = res[!is.na(res$padj),]
res = res[order(res$padj),]
print(paste(length(res[res$padj < 0.05 & res$log2FoldChange > 1,]$gene), "upregulated genes"))

upregulated_genes <- res[res$padj < 0.05 & res$log2FoldChange > 1,]$gene
GO_results <- enrichGO(gene = upregulated_genes,OrgDb = "org.Hs.eg.db",keyType = "ENSEMBL",ont = "ALL")
GO_results <- as.data.frame(GO_results)

top_GO_results <- GO_results %>%
  slice_head(n = 20)  # Change 'n' to the desired number of top terms

ggplot(top_GO_results, aes(x = reorder(Description, -p.adjust), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red", name = "Adjusted p-value") +
  labs(x = "", y = "", title = "GO Term Enrichment: Doxycycline-Induced Upregulation in MCF10A MYC Cells") +
  theme_minimal()

# MYC vs TM overexspression

# MYC vs TM overexpression:
coldata = coldata[c(2,4,6,8,10,12),]
counts = counts[,c(2,4,6,8,10,12)]
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~gene)

dds <- DESeq(dds)
resultsNames(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("gene", "conc"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p1 <- ggplot(pcaData, aes(PC1, PC2, color=conc, shape=gene)) +
  geom_point(size=3) +
  geom_text_repel(aes(label=name)) +  # Add non-overlapping labels from 'name'
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+
  theme_bw(base_size = 15)
p1

res = results(dds, contrast=c("gene","TM","MYC"))
res = as.data.frame(res)
res$gene = add_genesymbol(rownames(res))
res = as.data.frame(res)
res = res[!is.na(res$padj),]
res = res[order(res$padj),]
print(paste(length(res[res$padj < 0.05 & res$log2FoldChange > 1,]$gene), "upregulated genes"))
print(paste(length(res[res$padj < 0.05 & res$log2FoldChange < -1,]$gene), "downregulated genes"))
#res[res$gene %in% c(“IL1A”,“IL1B”,“IL1RN”,“AREG”,“CD74”,“GADD45B”,“BHLHE40”,“KDM6B”,“NFKB1”,“JUN”,“B2M”,“CASP4",“NLRP3”,“HLA-DRB1"),]
options(repr.plot.width=6,repr.plot.height=8,repr.plot.res = 300)
#Visualise ageing and trained immunity genes
de_mat = assay(vsd)
rownames(de_mat)= add_genesymbol(rownames(de_mat))
de_mat = as.data.frame(de_mat)
de_mat = de_mat[rownames(de_mat) %in% c(res[res$padj < 0.05 & res$log2FoldChange > 1.2  | res$log2FoldChange < -1.2,]$gene),]
anno = coldata[,c(1,2)]

# Volcano Plot
res$label <- ifelse((res$log2FoldChange > 1.5 | res$log2FoldChange < -1.5) & res$padj < 0.05, res$gene, NA)
p_threshold <- -log10(0.05)



ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_text_repel(aes(label = label), max.overlaps = Inf) + # Add labels for significant genes
  geom_hline(yintercept = p_threshold, linetype = "dashed", color = "red") +  # Dashed line for p < 0.05
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "blue") +   # Dashed lines for log2FoldChange thresholds
  theme_bw(base_size = 15)
ggsave("volcano_TM_0_1000_log2FC_1.pdf",width = 12, height = 9)
