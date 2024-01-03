library(edgeR)


setwd("Htseq_Count/")
# Function to perform edgeR analysis
perform_edgeR_analysis <- function(subset_counts, conditions, postfix) {
  # Create DGEList object
  y <- DGEList(counts = subset_counts, group = conditions)
  
  # Normalize counts
  y <- calcNormFactors(y)
  
  # Estimate common dispersion
  y <- estimateGLMCommonDisp(y, verbose = TRUE)
  cat("Common dispersion estimate:", y$common.dispersion, "\n")
  
  # Estimate trended dispersion (without specifying trend argument)
  y <- estimateGLMTrendedDisp(y)
  cat("Trended dispersion estimate:", y$trended.dispersion, "\n")
  
  # Estimate tagwise dispersion (without verbose argument)
  y <- estimateGLMTagwiseDisp(y)
  cat("Tagwise dispersions estimates (first 10):", head(y$tagwise.dispersion, 10), "\n")
  
  # Perform exact test
  et <- exactTest(y)
  
  # Extract p-values
  p_values <- et$table$PValue
  
  # Adjust p-values using Benjamini-Hochberg procedure
  padj_values_fdr <- p.adjust(p_values, method = "fdr")
  
  # Check the new padj_values_fdr
  summary(padj_values_fdr)
  
  # Add adjusted p-values to the exact test results table
  et$table$padj <- padj_values_fdr  # Corrected assignment
  
  # Extract top differentially expressed genes with adjusted p-values
  # top_tags <- head(et$table[order(et$table$padj), ], n = 10)
  
  # Process and save results as a data frame
  results <- as.data.frame(et$table)
  results$gene_name <- rownames(results)
  rownames(results) <- NULL
  assign(paste0("et_", postfix, "_df"), results, envir = .GlobalEnv)
  
  return(list(p_values = p_values, padj_values = padj_values_fdr))
}



# Data processing
conditions_551_552 <- metadata_C551_C552$Result
conditions_557_558 <- metadata_C557_C558$Result
conditions_557_559 <- metadata_C557_C559$Result
conditions_558_559 <- metadata_C558_C559$Result


# edgeR analysis for each condition
results_551_552 <- perform_edgeR_analysis(subset_551_552, conditions_551_552, "551_552")
results_557_558 <- perform_edgeR_analysis(subset_557_558, conditions_557_558, "557_558")
results_557_559 <- perform_edgeR_analysis(subset_557_559, conditions_557_559, "557_559")
results_558_559 <- perform_edgeR_analysis(subset_558_559, conditions_558_559, "558_559")


# Define FPKM and FC thresholds
fpkm_threshold <- 1
fc_threshold <- 2

# Filter based on FPKM and FC
filtered_results_551_552 <- et_551_552_df[et_551_552_df$Mean_C551_C552 > fpkm_threshold & abs(et_551_552_df$logFC) > fc_threshold, ]
filtered_results_557_558 <- et_557_558_df[et_557_558_df$Mean_C557_C558 > fpkm_threshold & abs(et_557_558_df$logFC) > fc_threshold, ]
filtered_results_557_559 <- et_557_559_df[et_557_559_df$Mean_C557_C559 > fpkm_threshold & abs(et_557_559_df$logFC) > fc_threshold, ]
filtered_results_558_559 <- et_558_559_df[et_558_559_df$Mean_C558_C559 > fpkm_threshold & abs(et_558_559_df$logFC) > fc_threshold, ]

# Save or further process filtered results as needed
significant_genes_551_552<-filtered_results_551_552$gene_name
significant_genes_557_558<-filtered_results_557_558$gene_name
significant_genes_557_559<-filtered_results_557_559$gene_name
significant_genes_558_559<-filtered_results_558_559$gene_name


# Create a list of significant genes at each time point
significant_genes_list <- list(S_551_552 = significant_genes_551_552,
                               S_557_558 = significant_genes_557_558,
                               S_557_559 = significant_genes_557_559,
                               S_558_559 = significant_genes_558_559 )


# Install and load the openxlsx package
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)
# Find the maximum length among all significant genes lists
max_length <- max(sapply(significant_genes_list, length))

# Pad each list with NA values to make them the same length
padded_genes_list <- lapply(significant_genes_list, function(x) c(x, rep(NA, max_length - length(x))))

# Create a data frame with columns named after each time point
significant_genes_df <- data.frame(padded_genes_list)

# Specify the file path where you want to save the Excel file
excel_file_path <- "FC2_significant_genes_dataframe.xlsx"

# Use write.xlsx to save the data frame to an Excel file
write.xlsx(significant_genes_df, file = excel_file_path, rowNames = FALSE)

