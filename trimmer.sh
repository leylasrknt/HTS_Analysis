#!/bin/bash

# Specify the adapter fasta file
adapter_fasta="path/to/your/adapter_file.fa"

# Specify the directory where your input files are located
INPUT_DIR="path/to/your/input/files"

# Specify the directory where your output files should be located
OUTPUT_DIR="path/to/your/output/files"

# Specify the directory where log files should be located
LOG_DIR="path/to/your/log/files"

# List of input file prefixes
FILE_PREFIXES=(
  "sample_1"
  "sample_2"
  "sample_3"
  "sample_4"
  "sample_5"
  "sample_6"
  "sample_7"
  "sample_8"
  "sample_9"
  "sample_10"
  ...
)

# Specify the Trimmomatic command (assuming it's in your PATH)
TRIMMER="ILLUMINACLIP:${adapter_fasta}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

for FP in "${FILE_PREFIXES[@]}"; do
  INFILES="$INPUT_DIR/${FP}_1.fastq $INPUT_DIR/${FP}_2.fastq"
  OUTFILES="$OUTPUT_DIR/${FP}_1_trimmed.fastq $OUTPUT_DIR/${FP}_1_unpaired.fastq $OUTPUT_DIR/${FP}_2_trimmed.fastq $OUTPUT_DIR/${FP}_2_unpaired.fastq"
  LOG="$LOG_DIR/trimmomatic/$FP.log"

  mkdir -p "$(dirname "$LOG")" "$OUTPUT_DIR"

  echo "Processing $FP..."
  echo "INFILES: $INFILES"
  echo "OUTFILES: $OUTFILES"
  echo "LOG: $LOG"

  trimmomatic PE -threads 2 $INFILES $OUTFILES $TRIMMER > "$LOG" 2>&1
done
