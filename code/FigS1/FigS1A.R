# Supplemental figure RNAseq (Fig S1A)
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
library(sva)

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
# First library prep for TP53
counts_1 = read.table("data/transcriptomics/counts/experiment_1/GRCh38.p14-counts.tsv",header = TRUE)
rownames(counts_1) = counts_1$gene
counts_1 <- counts_1 %>%
  dplyr::select(c(1:4,8:10,11:13,17:19,20:22,26:28))

# Second library prep for other conditions:
counts_2 = read.table("data/transcriptomics/counts/experiment_2/GRCh38.p14-counts.tsv",header = TRUE)
rownames(counts_2) = counts_2$gene

counts <- merge(counts_1, counts_2, by = "gene")
colnames(counts) <- c("gene","TP53R175H_0_1","TP53R175H_0_2","TP53R175H_0_3", "TP53R175H_1000_1","TP53R175H_1000_2","TP53R175H_1000_3","MYCbatch1_0_1","MYCbatch1_0_2","MYCbatch1_0_3","MYCbatch1_1000_1","MYCbatch1_1000_2","MYCbatch1_1000_3","EV_0_1","EV_0_2","EV_0_3","EV_1000_1","EV_1000_2", "EV_1000_3","TM_0_1","TM_1000_1","MYC_0_1","MYC_1000_1","TM_0_2","TM_1000_2","MYC_0_2","MYC_1000_2","TM_0_3","TM_1000_3","MYC_0_3","MYC_1000_3")
counts$gene = NULL
listmeta = as.data.frame(colnames(counts))

#Generate metadata file
df = as.data.frame(strsplit(listmeta[,1],split  = "_")[[1]])
df = as.data.frame(t(df))
for(i in 1:30){
  if(i == 1){
    next
  }
  suppressWarnings({
    df = rbind(df,strsplit(listmeta[,1],split  = "_")[[i]])
  })
}
df$batch <- c(rep("batch1", 18), rep("batch2", 12))
df$batch <- factor(df$batch)


colnames(df) = c("gene","conc","rep","batch")
coldata = df
coldata
rownames(coldata) = listmeta[,1]

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~batch +conc)
dds <- DESeq(dds)
resultsNames(dds)
vsd <- vst(dds, blind=FALSE)
vsd_bc <- limma::removeBatchEffect(assay(vsd), batch = coldata$batch, design = model.matrix(~ conc, data=coldata))
#rld <- rlog(dds, blind=FALSE)
#pcaData <- plotPCA(rld, intgroup=c("gene", "conc"), returnData=TRUE)
vsd2 <- vsd
assay(vsd2) <- vsd_bc
vsd_uncorrected <- vsd

# PCA
pcaData <- plotPCA(vsd_uncorrected, intgroup=c("gene", "conc"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$sample <- rownames(coldata)
pcaData$batch <- coldata$batch


p1 <- ggplot(pcaData, aes(PC1, PC2, color=batch, shape=conc)) +
  geom_point(size=2) +
  geom_text_repel(aes(label=sample), size=5, max.overlaps = 20) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw(base_size =14)
p1
#ggsave(filename = "../figures/final_figures/Fig2/supplemental/pca_plot_RNAseq_allsamples_uncorrected.pdf", plot = p1, width = 8, height = 6)

# With combat:
batch <- coldata$batch
mod <- model.matrix(~ conc, data=coldata)
vsd_mat <- assay(vsd)
vsd_combat <- ComBat(dat=vsd_mat, batch=batch, mod=mod)


# Remove TM and MYC conditions
vsd_combat <-vsd_combat %>%as.data.frame() %>%
  dplyr::select(
    TP53R175H_0_1, TP53R175H_0_2, TP53R175H_0_3,
    TP53R175H_1000_1, TP53R175H_1000_2, TP53R175H_1000_3,
    EV_0_1, EV_0_2, EV_0_3,
    EV_1000_1, EV_1000_2, EV_1000_3
  )
# Run PCA
pca_res <- prcomp(t(vsd_combat))
coldata <- coldata %>% filter(gene == "TP53R175H" | gene == "EV")

# Extract PC1 and PC2 from prcomp result
pca_data <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  gene = coldata$gene,
  conc = coldata$conc,
  batch = coldata$batch
)
pca_data$sample <- rownames(coldata)

# Calculate percent variance explained
percentVar <- round(100 * (pca_res$sdev^2) / sum(pca_res$sdev^2))

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = conc, shape = gene)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = sample), size = 6, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 15)

# Remove sample MYCmut_0_3 from pca_data
pca_data_filtered <- subset(pca_data, gene != c("MYCbatch1"))

p2<- ggplot(pca_data_filtered, aes(x = PC1, y = PC2, color = conc, shape = gene)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = sample), size = 5, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14)
p2
#ggsave(filename = "../figures/final_figures/Fig2/supplemental/pca_plot_RNAseq_allsamples_batchcorrected.pdf", plot = p2, width = 8, height = 6)












