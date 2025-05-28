#!/bin/bash

# Parse command-line arguments
INPUT_DIR=${1:-"./01_data"}
OUTPUT_DIR=${2:-"./02_QC/fastqc"}
REPORT_DIR=${3:-"./02_QC/_multiqc"}

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$REPORT_DIR"

# Run fastQC on all paired-end reads
for file in "$INPUT_DIR"/*_R1.fastq.gz; do
    sample=$(basename "$file" _R1.fastq.gz)
    echo "Running fastQC for $sample"
    fastqc "$INPUT_DIR/${sample}_R1.fastq.gz" "$INPUT_DIR/${sample}_R2.fastq.gz" -o "$OUTPUT_DIR"
done

# Run multiQC on all FastQC reports
echo "Running multiQC"
multiqc "$OUTPUT_DIR" -o "$REPORT_DIR"

echo "QC completed. Reports are in $REPORT_DIR"
