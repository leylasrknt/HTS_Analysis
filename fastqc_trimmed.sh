#!/bin/bash

# Activate the conda environment named "fastqc"
conda activate fastqc

# List of trimmed fastq files (replace these with your actual filenames)
fastq_files=(
    "sample1_trimmed.fastq"
    "sample2_trimmed.fastq"
    "sample3_trimmed.fastq"
    "sample4_trimmed.fastq"
)

# Output directory for FastQC results (replace this with your desired output path)
output_directory="/path/to/output_directory"

# Run FastQC for each trimmed fastq file
for file in "${fastq_files[@]}"; do
    # Run FastQC on the current file and save the results to the specified output directory
    fastqc "$file" -o "$output_directory"
done
