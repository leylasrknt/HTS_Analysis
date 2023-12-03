#!/bin/bash

# Directory containing your BAM files
bam_dir="/path/to/your/bam/directory"

# Directory to save the Picard results
output_dir="/path/to/your/output/directory/Picard_Results_Optical_Duplication"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each BAM file in the specified directory
for bam_file in "$bam_dir"/*.sorted.bam; do
    # Check if there are any matching files
    if [ -e "$bam_file" ]; then
        # Define the output file name (with the output directory)
        output_bam="$output_dir/$(basename "${bam_file%.sorted.bam}.marked.bam")"
        # Define the metrics file name (with the output directory)
        metrics_file="$output_dir/$(basename "${bam_file%.sorted.bam}.metrics.txt")"

        echo "Processing $bam_file..."

        # Run Picard's MarkDuplicates command
        picard MarkDuplicates \
            INPUT="$bam_file" \
            OUTPUT="$output_bam" \
            METRICS_FILE="$metrics_file" \
            REMOVE_DUPLICATES=true

        echo "Done processing $bam_file"
    else
        echo "No BAM files found in $bam_dir"
    fi
done
