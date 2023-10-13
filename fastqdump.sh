#!/bin/bash

# Path to fasterq-dump tool
fastqDump="/bioinf/tools/sratoolkit.3.0.6-ubuntu64/bin/fasterq-dump"

# Command to download and split SRA file into fastq files
$sraFile="your_sra_file.sra" # Replace with the actual SRA file name
$fastqDir="/path/to/your/output_directory" # Replace with the desired output directory

# Run fasterq-dump to download and split SRA file
$fastqDump --split-3 --outdir $fastqDir $sraFile
