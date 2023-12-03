#!/bin/bash

# Set the paths and variables
INDEX_DIR="path/to/your/genome_index"
OUTPUT_DIR="path/to/your/output/directory"
FASTQ_DIR="path/to/your/trimmed/fastq/files"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# List of sample numbers
SAMPLES=("abbreviations")

# Get the list of trimmed FASTQ files in the specified directory
# Loop through each sample and perform alignment
for i in "${SAMPLES[@]}"; do
    echo "Processing SRR99375$i..."

    # Define the file names for the current sample
    TRIMMED_R1="${FASTQ_DIR}/SRR99375${i}_1_trimmed.fastq"
    TRIMMED_R2="${FASTQ_DIR}/SRR99375${i}_2_trimmed.fastq"

    # Run HISAT2 for paired-end data
    hisat2 -x "$INDEX_DIR/genome_index" -1 "$TRIMMED_R1" -2 "$TRIMMED_R2" -S "$OUTPUT_DIR/SRR99375${i}.sam"

    # Convert SAM to BAM
    samtools view -b -o "$OUTPUT_DIR/SRR99375${i}.bam" "$OUTPUT_DIR/SRR99375${i}.sam"

    # Sort and index the BAM file
    samtools sort -o "$OUTPUT_DIR/SRR99375${i}.sorted.bam" "$OUTPUT_DIR/SRR99375${i}.bam"
    samtools index "$OUTPUT_DIR/SRR99375${i}.sorted.bam"

    # Remove intermediate files (SAM and unsorted BAM)
    rm "$OUTPUT_DIR/SRR99375${i}.sam"
    rm "$OUTPUT_DIR/SRR99375${i}.bam"

    echo "SRR99375$i completed."
done

echo "All alignment jobs completed."
