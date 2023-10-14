# Load required libraries
library(clusterProfiler)
library(biomaRt)
library(dplyr)

# Specify the dataset and organism you are working with
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to convert gene names to ENSEMBL IDs and perform enrichment analysis
perform_enrichment_analysis <- function(gene_list) {
  # Convert gene names to ENSEMBL IDs
  gene_ids <- getBM(attributes = c("ensembl_gene_id"), 
                    filters = "external_gene_name", 
                    values = gene_list, 
                    mart = ensembl)
  
  # Extract the ENSEMBL IDs from the result
  your_gene_list <- gene_ids$ensembl_gene_id
  
  # Perform gene set enrichment analysis
  go_results <- enrichGO(gene = your_gene_list, 
                         keyType = "ENSEMBL", 
                         OrgDb = org.Hs.eg.db,
                         ont = "BP", 
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2)
  
  return(go_results)
}

# Gene lists for upregulated genes (replace these with your actual data)
gene_lists_upregulated <- list(
  upregulated_genes_4,
  upregulated_genes_8,
  upregulated_genes_12,
  upregulated_genes_24,
  upregulated_genes_48,
  upregulated_genes_168
)

# Gene lists for downregulated genes (replace these with your actual data)
gene_lists_downregulated <- list(
  downregulated_genes_4,
  downregulated_genes_8,
  downregulated_genes_12,
  downregulated_genes_24,
  downregulated_genes_48,
  downregulated_genes_168
)

# Perform enrichment analysis for each gene list in gene_lists_upregulated
enrichment_results_upregulated <- lapply(gene_lists_upregulated, perform_enrichment_analysis)

# Perform enrichment analysis for each gene list in gene_lists_downregulated
enrichment_results_downregulated <- lapply(gene_lists_downregulated, perform_enrichment_analysis)
