#!/bin/bash

# Set the paths and variables
INDEX_DIR="path/to/your/genome_index"
OUTPUT_DIR="path/to/your/output/directory"
FASTQ_DIR="path/to/your/trimmed/fastq/files"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Get the list of trimmed FASTQ files in the specified directory
TRIMMED_FILES=("$FASTQ_DIR"/*_1.fastq.gz)

# Loop through each pair of trimmed FASTQ files and perform alignment
for trimmed_file_1 in "${TRIMMED_FILES[@]}"; do
    # Assuming the files are named consistently with "_1" and "_2" for the first and second reads
    trimmed_file_2="${trimmed_file_1/_1/_2}"

    echo "Processing $trimmed_file_1 and $trimmed_file_2..."

    # Extract the sample identifier from the file name (modify as needed)
    sample_id=$(basename "$trimmed_file_1" | sed 's/sample_\([0-9]*\)_1.*\.fastq/\1/')

    # Run HISAT2 for paired-end data
    hisat2 -x "$INDEX_DIR/genome_index" -1 "$trimmed_file_1" -2 "$trimmed_file_2" -S "$OUTPUT_DIR/sample_${sample_id}.sam"

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
