# Install and load UpSetR package
install.packages("UpSetR")
library(UpSetR)

# Function to create gene data matrix and generate upset plot
create_upset_plot <- function(gene_lists, color) {
  # Find the union of all genes
  all_genes <- unique(unlist(gene_lists))
  
  # Create an empty matrix
  gene_matrix <- matrix(0, nrow = length(all_genes), ncol = length(gene_lists))
  
  # Fill in the matrix with 1 for present genes
  for (i in 1:length(gene_lists)) {
    present_genes <- gene_lists[[i]]
    gene_matrix[all_genes %in% present_genes, i] <- 1
  }
  
  # Create a data frame from the matrix
  gene_data <- as.data.frame(gene_matrix)
  colnames(gene_data) <- paste("Set", seq_along(gene_lists), sep = "")
  
  # Create the upset plot
  upset(gene_data, sets.bar.color = color, order.by = "freq")
}

# Combine all gene vectors into a list for upregulated genes
gene_lists_upregulated <- list(
  upregulated_genes_4,
  upregulated_genes_8,
  upregulated_genes_12,
  upregulated_genes_24,
  upregulated_genes_48,
  upregulated_genes_168
)

# Generate upset plot for upregulated genes
create_upset_plot(gene_lists_upregulated, "#56B4E9")

# Combine all gene vectors into a list for downregulated genes
gene_lists_downregulated <- list(
  downregulated_genes_4,
  downregulated_genes_8,
  downregulated_genes_12,
  downregulated_genes_24,
  downregulated_genes_48,
  downregulated_genes_168
)

# Generate upset plot for downregulated genes
create_upset_plot(gene_lists_downregulated, "#D55E00")
