#!/bin/bash

# Directory containing your BAM files
BAM_DIR="path/to/your/BAM/files"

# Ensure the output directory exists
INDEX_DIR="path/to/your/index/directory"
mkdir -p "$INDEX_DIR"

# Loop through each BAM file in the directory and create index (BAI) files
for BAM_FILE in "$BAM_DIR"/*.bam; do
  if [ -f "$BAM_FILE" ]; then
    # Check if the file is a regular file
    echo "Indexing $BAM_FILE"
    # Create output file name for BAI
    INDEX_FILE="$INDEX_DIR/$(basename "$BAM_FILE").bai"
    samtools index "$BAM_FILE" "$INDEX_FILE"
  fi
done
