# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tidyverse)

# Create a data frame for plotting
plot_data <- data.frame(
  log2FoldChange = results$log2FoldChange,
  padj = results$padj,
  Significant = ifelse(results$padj < 0.05, "Significant", "Not Significant")
)

# Volcano plot
volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point(size = 2) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 (Adjusted p-value)") +
  theme_minimal()

# PCA plot
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "TimePoint"), returnData=TRUE)
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Treatment, shape = TimePoint)) +
  geom_point(size = 3) +
  ggtitle("PCA Analysis of DESeq Data") +
  theme_minimal()

# Heatmap
vsd2 <- na.omit(vsd)
topVarGenes <- head(order(rowVars(assay(vsd2)), decreasing = TRUE), 50)
topVarGenes1 <- as.data.frame(assay(vsd)[topVarGenes, ])
pheatmap(topVarGenes1, scale = "row",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

# Filter up-regulated and down-regulated genes for each time point
filter_genes <- function(results_df, threshold, padj_threshold) {
  filtered_results <- na.omit(results_df)
  upregulated <- filtered_results[filtered_results$padj < padj_threshold & filtered_results$log2FoldChange > threshold, ]
  downregulated <- filtered_results[filtered_results$padj < padj_threshold & filtered_results$log2FoldChange < -threshold, ]
  return(list(upregulated, downregulated))
}

up_down_results <- function(results_object, padj_threshold = 0.05) {
  results_df <- as.data.frame(results(results_object))
  list(upregulated = filter_genes(results_df, 1, padj_threshold),
       downregulated = filter_genes(results_df, 1, padj_threshold))
}

# Example usage for time points
up_down_results_4 <- up_down_results(dds_4)
upregulated_genes_4 <- up_down_results_4$upregulated$gene_name
downregulated_genes_4 <- up_down_results_4$downregulated$gene_name
