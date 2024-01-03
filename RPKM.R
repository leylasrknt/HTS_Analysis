# Load necessary libraries
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library("GenomicFeatures")
library(Rsubread)
library(dplyr)

# Path to your GTF file
gtf_file <- "/genomic.gtf"

# Read GTF file
gtf <- readGFF(gtf_file)###subset 

###############################################################################
exon_gtf <- gtf[gtf$type == "exon", ]
exon_lengths <- tapply(exon_gtf$end - exon_gtf$start + 1, exon_gtf$gene_id, sum)
exon_lengths_df <- data.frame(Gene = names(exon_lengths), Length = as.numeric(exon_lengths))
###############################################################################
gene_length <- subset(gtf, type == "exon", select = c(gene, start, end))
gene_length 

# Calculate gene lengths
gene_length$length <- (gene_length$end - gene_length$start) + 1


#######


# Calculate RPKM values for each gene
# Extract gene names from combined_df
gene_names <- rownames(combined_df)
col<-colnames(combined_df)
# Filter gene_length dataframe to include only rows with matching gene names
filtered_gene_length <- exon_lengths_df[exon_lengths_df$Gene %in% gene_names, ]

# Initialize a data frame to store RPKM values
rpkm_df <- data.frame(matrix(NA, nrow = nrow(combined_df), ncol = ncol(combined_df)))

# Set row names of rpkm_df to gene names from combined_df
rownames(rpkm_df) <- gene_names
colnames(rpkm_df) <- col

# Calculate RPKM values for each column and store them in the new data frame
# Calculate RPKM values for each gene and store them in rpkm_df
for (gene_name in gene_names) {
  for (sample_col in col) {
    # Calculate RPKM for the current gene and sample
    reads <- combined_df[gene_name, sample_col]
    exon_length <- filtered_gene_length[filtered_gene_length$Gene == gene_name, "Length"]
    
    rpkm <- (reads * 1e9) / (exon_length * sum(reads))
    
    # Assign the calculated RPKM value to rpkm_df
    rpkm_df[gene_name, sample_col] <- rpkm
  }
}

# View rpkm_df
print(rpkm_df)
###############################################################################
#second way
# Calculate RPKM values for each gene and store them in rpkm_df
# Calculate RPKM values for each gene and store them in the rpkm_df data frame
for (gene_name in gene_names) {
  gene_counts <- combined_df[gene_name, ]
  gene_length <- filtered_gene_length$Length[filtered_gene_length$Gene == gene_name]
  
  # Calculate RPKM for each sample
  rpkm_values <- (gene_counts * 1e9) / (gene_length * sum(gene_counts))
  
  # Assign RPKM values to the corresponding row in rpkm_df
  rpkm_df[gene_name, ] <- rpkm_values
}


# View rpkm_df
print(rpkm_df)
write.table( rpkm_df,"genomic.gtf/rpkm.csv", sep = ",", col.names = TRUE, row.names = TRUE)



