#!/bin/bash

# Set the paths and variables
INDEX_DIR="path/to/your/genome_index"
OUTPUT_DIR="path/to/your/output/directory"
FASTQ_DIR="path/to/your/trimmed/fastq/files"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# List of sample numbers or identifiers
SAMPLES=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24")

# Loop through each sample and perform alignment
for i in "${SAMPLES[@]}"; do
    echo "Processing Sample $i..."

    # Define the file names for the current sample
    TRIMMED_R1="${FASTQ_DIR}/sample_${i}_1_trimmed.fastq"
    TRIMMED_R2="${FASTQ_DIR}/sample_${i}_2_trimmed.fastq"

    # Run HISAT2 for paired-end data
    hisat2 -x "$INDEX_DIR/genome_index" -1 "$TRIMMED_R1" -2 "$TRIMMED_R2" -S "$OUTPUT_DIR/sample_${i}.sam"

    # Convert SAM to BAM
    samtools view -b -o "$OUTPUT_DIR/sample_${i}.bam" "$OUTPUT_DIR/sample_${i}.sam"

    # Sort and index the BAM file
    samtools sort -o "$OUTPUT_DIR/sample_${i}.sorted.bam" "$OUTPUT_DIR/sample_${i}.bam"
    samtools index "$OUTPUT_DIR/sample_${i}.sorted.bam"

    # Remove intermediate files (SAM and unsorted BAM)
    rm "$OUTPUT_DIR/sample_${i}.sam"
    rm "$OUTPUT_DIR/sample_${i}.bam"

    echo "Sample $i completed."
done

echo "All alignment jobs completed."
