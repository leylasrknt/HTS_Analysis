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


########################
# Load necessary libraries
library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library("GenomicFeatures")
library(Rsubread)
library(dplyr)

# Path to your GTF file
gtf_file <- "genomic.gtf"

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
write.table( rpkm_df,"rpkm.csv", sep = ",", col.names = TRUE, row.names = TRUE)
###############################################################################
for (gene_name in gene_names) {
  gene_counts <- combined_df[gene_name, ]
  gene_length <- filtered_gene_length$Length[filtered_gene_length$Gene == gene_name]
  
  cat("Gene Name:", gene_name, "\n")
  print(gene_counts)
  cat("Gene Length:", gene_length, "\n")
  
  # Calculate RPKM for each sample
  denominator <- sum(as.numeric(gene_counts))
  
  # Check if the denominator is zero
  if (denominator == 0) {
    rpkm_values <- numeric(length(gene_counts))
  } else {
    rpkm_values <- as.numeric(gene_counts) * 1e9 / (as.numeric(gene_length) * denominator)
  }
  
  cat("RPKM Values:")
  print(rpkm_values)
  cat("\n")
  
  # Assign RPKM values to the corresponding row in rpkm_df
  if (length(rpkm_values) == 0) {
    # Handle the case where rpkm_values is of length 0
    rpkm_df[gene_name, ] <- numeric(length(gene_counts))
  } else {
    rpkm_df[gene_name, ] <- rpkm_values
  }
}

# Print rpkm_df
print(rpkm_df)
###############################################################################

###############################################################################
# Subset the rpkm_df for the specific samples and matching gene names
rpkm_subset_C551_S3_C552_S4 <- subset(rpkm_df, select = c(C551_S3, C552_S4))
rpkm_subset_C551_S3_C552_S4 <- rpkm_subset_C551_S3_C552_S4[intersect(rownames(rpkm_df), et_551_552_df$gene_name),]

# Add the RPKM values as new columns in results
et_551_552_df$C551_S3 <- rpkm_subset_C551_S3_C552_S4$C551_S3
et_551_552_df$C552_S4 <- rpkm_subset_C551_S3_C552_S4$C552_S4


##############################################################################

# Subset the rpkm_df for the specific samples and matching gene names
rpkm_subset_C557_S5_C558_S6 <- subset(rpkm_df, select = c(C557_S5,C558_S6))
rpkm_subset_C557_S5_C558_S6 <- rpkm_subset_C557_S5_C558_S6[intersect(rownames(rpkm_df), et_557_558_df$gene_name),]

# Add the RPKM values as new columns in results
et_557_558_df$C557_S5 <- rpkm_subset_C557_S5_C558_S6$C557_S5
et_557_558_df$C558_S6 <- rpkm_subset_C557_S5_C558_S6$C558_S6


##############################################################################

# Subset the rpkm_df for the specific samples and matching gene names
rpkm_subset_C557_S5_C559_S7 <- subset(rpkm_df, select = c(C557_S5,C559_S7))
rpkm_subset_C557_S5_C559_S7 <- rpkm_subset_C557_S5_C559_S7[intersect(rownames(rpkm_df), et_557_559_df$gene_name),]

# Add the RPKM values as new columns in results
et_557_559_df$C557_S5 <- rpkm_subset_C557_S5_C559_S7$C557_S5
et_557_559_df$C559_S7 <- rpkm_subset_C557_S5_C559_S7$C559_S7



##############################################################################

# Subset the rpkm_df for the specific samples and matching gene names
rpkm_subset_C558_S6_C559_S7 <- subset(rpkm_df, select = c(C558_S6,C559_S7))
rpkm_subset_C558_S6_C559_S7 <- rpkm_subset_C558_S6_C559_S7[intersect(rownames(rpkm_df), et_558_559_df$gene_name),]

# Add the RPKM values as new columns in results
et_558_559_df$C558_S6 <- rpkm_subset_C558_S6_C559_S7$C558_S6
et_558_559_df$C559_S7 <- rpkm_subset_C558_S6_C559_S7$C559_S7


# Calculate the mean of C551_S3 and C552_S4 for et_551_552_df
et_551_552_df$Mean_C551_C552 <- rowMeans(rpkm_subset_C551_S3_C552_S4[, c("C551_S3", "C552_S4")], na.rm = TRUE)

# Calculate the mean of C557_S5 and C558_S6 for et_557_558_df
et_557_558_df$Mean_C557_C558 <- rowMeans(rpkm_subset_C557_S5_C558_S6[, c("C557_S5", "C558_S6")], na.rm = TRUE)

# Calculate the mean of C557_S5 and C559_S7 for et_557_559_df
et_557_559_df$Mean_C557_C559 <- rowMeans(rpkm_subset_C557_S5_C559_S7[, c("C557_S5", "C559_S7")], na.rm = TRUE)

# Calculate the mean of C558_S6 and C559_S7 for et_558_559_df
et_558_559_df$Mean_C558_C559 <- rowMeans(rpkm_subset_C558_S6_C559_S7[, c("C558_S6", "C559_S7")], na.rm = TRUE)



