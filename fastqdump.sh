#!/bin/bash

# Path to fasterq-dump tool
fastqDump="/pathwaytofastqdump"

# Command to download and split SRA file into fastq files
$sraFile="your_sra_file.sra" # Replace with the actual SRA file name
$fastqDir="/path/to/your/output_directory" # Replace with the desired output directory

# Run fasterq-dump to download and split SRA file
$fastqDump --split-3 --outdir $fastqDir $sraFile

##loop for multiple SRA file
# Output directory
outputDir="/path/to/your/output_directory" # Replace with your desired output directory

# Read SRA accession numbers from the file and download FASTQ files
while IFS= read -r accession; do
    $fastqDump --split-3 --outdir "$outputDir" "$accession"
done < sra_accessions.txt

echo "Download of all SRA files completed."
