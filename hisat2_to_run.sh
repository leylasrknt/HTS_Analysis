#!/bin/bash

# Set the paths and variables
INDEX_DIR="path/to/your/genome_index"
OUTPUT_DIR="path/to/your/output/directory"
FASTQ_DIR="path/to/your/trimmed/fastq/files"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Get the list of trimmed FASTQ files in the specified directory
TRIMMED_FILES=("$FASTQ_DIR"/*.fastq)

# Loop through each trimmed FASTQ file and perform alignment
for trimmed_file in "${TRIMMED_FILES[@]}"; do
    echo "Processing $trimmed_file..."

    # Extract the sample identifier from the file name (modify as needed)
    sample_id=$(basename "$trimmed_file" | sed 's/sample_\([0-9]*\)_.*\.fastq/\1/')

    # Run HISAT2 for paired-end data
    hisat2 -x "$INDEX_DIR/genome_index" -U "$trimmed_file" -S "$OUTPUT_DIR/sample_${sample_id}.sam"

    # Convert SAM to BAM
    samtools view -b -o "$OUTPUT_DIR/sample_${sample_id}.bam" "$OUTPUT_DIR/sample_${sample_id}.sam"

    # Sort and index the BAM file
    samtools sort -o "$OUTPUT_DIR/sample_${sample_id}.sorted.bam" "$OUTPUT_DIR/sample_${sample_id}.bam"
    samtools index "$OUTPUT_DIR/sample_${sample_id}.sorted.bam"

    # Remove intermediate files (SAM and unsorted BAM)
    rm "$OUTPUT_DIR/sample_${sample_id}.sam"
    rm "$OUTPUT_DIR/sample_${sample_id}.bam"

    echo "Sample $sample_id completed."
done

echo "All alignment jobs completed."
