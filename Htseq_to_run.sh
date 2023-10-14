#!/bin/bash

# Path to the annotation file (GFF/GTF)
ANNOTATION_FILE="/path/to/your/annotation_file.gff"

# Output directory for count files
OUTPUT_DIR="/path/to/your/output/directory"

# Directory containing your BAM files
BAM_DIR="/path/to/your/bam/files"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Create an array of BAM files by listing them in the directory
bam_files=("$BAM_DIR"/*.bam)

# Loop through each BAM file in the bam_files array and run htseq-count
for BAM_FILE in "${bam_files[@]}"; do
  # Extract the file name (without extension) for naming the output count file
  FILENAME=$(basename "$BAM_FILE" .bam)

  # Run htseq-count
  htseq-count -f bam -r pos -t exon -i gene_id "$BAM_FILE" "$ANNOTATION_FILE" > "$OUTPUT_DIR/${FILENAME}_counts.txt"

  # Check if htseq-count was successful
  if [ $? -eq 0 ]; then
    echo "HTSeq counting for ${FILENAME} completed successfully."
  else
    echo "Error: HTSeq counting for ${FILENAME} failed."
  fi
done

echo "HTSeq counting for all samples completed."
