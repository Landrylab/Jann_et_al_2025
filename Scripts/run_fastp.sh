#!/bin/bash
# Set input and output directories from command-line arguments
INPUT_DIR=${1:-"./01_data"}
OUTPUT_DIR=${2:-"./03_filter"}
REPORT_DIR="${OUTPUT_DIR}/reports"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR" "$REPORT_DIR"

# Set R1 file pattern
R1_PATTERN="*_R1.fastq.gz"

# Loop through R1 files
for read in ${INPUT_DIR}/${R1_PATTERN}; do
  echo "Processing: $read"
  
  # Extract sample name by removing path and suffix
  sample=$(basename "$read" | sed -E 's/_R1\.fastq\.gz//')

  # Define the corresponding R2 file
  R2="${INPUT_DIR}/${sample}_R2.fastq.gz"

  # Check if R2 exists before proceeding
  if [[ -f "$R2" ]]; then
    echo "Found R2 for $sample: $R2"
    
    fastp \
      --in1 "$read" \
      --in2 "$R2" \
      --qualified_quality_phred 15 \
      --unqualified_percent_limit 40 \
      --out1 "${OUTPUT_DIR}/${sample}_R1.fastq.gz" \
      --out2 "${OUTPUT_DIR}/${sample}_R2.fastq.gz" \
      --html "${REPORT_DIR}/${sample}_fastp.html" \
      --json "${REPORT_DIR}/${sample}_fastp.json"
  else
    echo "Warning: R2 file not found for $sample. Skipping..."
  fi
done

# Run MultiQC on fastp reports
echo "Running MultiQC on fastp reports..."
multiqc "$REPORT_DIR" -o "$REPORT_DIR"

echo "Quality control pipeline completed."

