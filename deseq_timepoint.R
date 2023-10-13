# Load required libraries
library(DESeq2)
library(dplyr)

# Set working directory
setwd("/path/to/your/count_data_directory")

# List all count files in the directory and store them in a list
count_files <- list.files(pattern = ".txt", full.names = TRUE)
count_list <- list()

# Loop through the count files, read them, and store in the count_list
for (count_file in count_files) {
  sample_name <- gsub("\\.txt", "", basename(count_file))
  count_data <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t")
  colnames(count_data) <- sample_name
  count_list[[sample_name]] <- count_data
}

# Read metadata (replace 'metadata.csv' with your metadata file name)
metadata <- read.csv("metadata.csv", header = TRUE)
metadata$Treatment <- as.factor(metadata$Treated)
metadata$TimePoint <- as.factor(metadata$TimePoint)

# Create a DESeqDataSet object with the updated metadata and design formula
dds <- DESeqDataSetFromMatrix(countData = do.call(cbind, count_list),
                              colData = metadata,
                              design = ~ Treatment + TimePoint)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get differential expression results
results <- results(dds)

# Save results as CSV files
write.csv(results, file = "differential_expression_results.csv")
