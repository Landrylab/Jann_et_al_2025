#!/bin/bash


#conda activate seqprep


#R1=../DATA/PDR1_GiF13_10_R1.fastq.gz
#R2=../DATA/PDR1_GiF13_10_R2.fastq.gz

#R1=./03_filter/Pool*_R1.fastq.gz

#for read in $R1
#  do
#  echo $read
#  fragment=$(echo $read | cut -d"/" -f4 | cut -d"_" -f2)
#  pcr=$(echo $read | cut -d"/" -f4 | cut -d"_" -f3)
#  bbmerge.sh in1=../05_filter/filter/PDR1_${fragment}_${pcr}_R1.fastq.gz in2=../05_filter/filter/PDR1_${fragment}_${pcr}_R2.fastq.gz out=./merged/PDR1_${fragment}_${pcr}.fasta outu1=./unmerged/PDR1_${fragment}_${pcr}_R1_unmerged.fasta outu2=./unmerged/PDR1_${fragment}_${pcr}_R2_unmerged.fasta -ordered=t
#  done

#for read in $R1
#  do
#  echo $read
#  sample=$(echo $read | cut -d"/" -f3 | cut -d"_" -f1)
#  bbmerge.sh in1=./03_filter/${sample}_R1.fastq.gz in2=./03_filter/${sample}_R2.fastq.gz out=./04_merged/${sample}.fasta outu1=./04_merged/unmerged/${sample}_R1_unmerged.fasta outu2=./04_merged/unmerged/${sample}_R2_unmerged.fasta -ordered=t
#  done

#conda deactivate
#bbmerge.sh in1=$R1 in2=$R2 out=./merged/test.fasta outu1=./merged/test_unmerged_r1.fasta outu2=./merged/test_unmerged_r2.fasta -ordered=t

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
