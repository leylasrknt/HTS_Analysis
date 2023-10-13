#!/bin/bash

# Directory containing your BAM files
bam_dir="/path/to/your/bam/directory"

# Directory to save the Picard results
output_dir="/path/to/your/output/directory/Picard_Results_Optical_Duplication"

# List of your BAM files (without the full path)
bam_files=(
    "file1.sorted.bam"
    "file2.sorted.bam"
    "file3.sorted.bam"
    # Add more BAM files as needed
)

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each BAM file and run Picard's MarkDuplicates command
for bam_file in "${bam_files[@]}"; do
    # Define the full path to the BAM file
    full_bam_path="$bam_dir/$bam_file"
    # Define the output file name (with the output directory)
    output_bam="$output_dir/$(basename "${full_bam_path%.sorted.bam}.marked.bam")"
    # Define the metrics file name (with the output directory)
    metrics_file="$output_dir/$(basename "${full_bam_path%.sorted.bam}.metrics.txt")"

    echo "Processing $full_bam_path..."

    # Run Picard's MarkDuplicates command
    picard MarkDuplicates \
        INPUT="$full_bam_path" \
        OUTPUT="$output_bam" \
        METRICS_FILE="$metrics_file" \
        REMOVE_DUPLICATES=true

    echo "Done processing $full_bam_path"
done
