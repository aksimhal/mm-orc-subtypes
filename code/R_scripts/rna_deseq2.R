# load library 
# DESeq2 version 1.34.0

library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)

library(tximport)
library(ggplot2)
library(ggrepel)

# # get count dataset
count_matrix <- read.csv("[insert file name]")
count_matrix = subset(count_matrix, select = -c(X) )


# # view first two rows
count_matrix <- round(count_matrix)

head(count_matrix)
samplelist = c(colnames(count_matrix) )

labellist <- read.csv("[insert file name]")

sampleTable <- data.frame(condition = factor(labellist$X0))
coldata <- data.frame(
  sample = samplelist,
  condition = sampleTable$condition, 
  row.names = "sample" )


coldata$condition <- as.factor(coldata$condition)

all(rownames(coldata) %in% colnames(count_matrix))
all(rownames(coldata) == colnames(count_matrix))


dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata,
                              design = ~ condition)

genes_to_remove = rowSums(counts(dds)) >= 10
write.csv(genes_to_remove, "genes_to_remove_rna_RNAtest.csv")


dds <- dds[rowSums(counts(dds)) >= 10,]

# # set control condition as reference
dds$condition <- relevel(dds$condition, ref = "upper")

dds <- DESeq(dds)

resultsNames(dds)

# # get gene expression table
# # at this step independent filtering is applied by default to remove low count genes
# # independent filtering can be turned off by passing independentFiltering=FALSE to results
res <- results(dds, alpha=0.01)  
# same as results(dds, name="condition_infected_vs_control") or 
# results(dds, contrast = c("condition", "infected", "control") )
summary(res) 

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot"))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 & abs(log2FoldChange)>3.5), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>3.5), points(log2FoldChange, -log10(padj), pch=20, col="red"))

log2FoldChange = res$log2FoldChange
padj = res$padj

save(log2FoldChange, padj, file = "rna_padj_log2FoldChange_RNAtest.RData")

