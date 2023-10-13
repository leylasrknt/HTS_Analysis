#!/bin/bash

# Define input directories (customize these paths according to your project structure)
fastqc_results="path/to/fastqc/results"
trimmomatic_logs="path/to/trimmomatic/logs"

# Output directory for MultiQC report (customize this path according to your desired output location)
multiqc_output="path/to/output/directory"

# Run MultiQC to generate a comprehensive QC report
multiqc -v -f -o "$multiqc_output" "$fastqc_results" "$trimmomatic_logs"
