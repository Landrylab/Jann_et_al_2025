#!/bin/bash
# Set input and output directories from command-line arguments
INPUT_DIR=${1:-"./03_filter"}
OUTPUT_DIR=${2:-"./04_merged"}
UNMERGED_DIR="${OUTPUT_DIR}/unmerged"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR" "$UNMERGED_DIR"

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
    
    bbmerge.sh \
      in1="$read" \
      in2="$R2" \
      out="${OUTPUT_DIR}/${sample}.fasta" \
      outu1="${UNMERGED_DIR}/${sample}_R1_unmerged.fasta" \
      outu2="${UNMERGED_DIR}/${sample}_R2_unmerged.fasta" \
      ordered=t
  else
    echo "Warning: R2 file not found for $sample. Skipping..."
  fi
done

echo "Merging pipeline completed."
