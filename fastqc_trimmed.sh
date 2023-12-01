#!/bin/bash

# Activate the conda environment named "fastqc"
conda activate fastqc

# Pattern for fastq files
fastq_pattern="/path/to/trimmed_files/*.fastq"

# Output directory for FastQC results 
output_directory="/path/to/output_directory"

# Run FastQC for each  fastq file matching the pattern
for file in $fastq_pattern; do
    # Run FastQC on the current file and save the results to the specified output directory
    fastqc "$file" -o "$output_directory"
done
