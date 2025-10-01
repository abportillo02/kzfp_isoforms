#!/bin/bash

# Check if input file with sample names is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <sample_list.txt>"
  exit 1
fi

samples="$1"
indir=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/Star_files
outdir=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/filtered_bams
samtools=/home/abportillo/.conda/envs/mamba_abner_BC/bin/samtools

# Create output directory if it doesn't exist
mkdir -p "$outdir"

while IFS= read -r sample_name; do
  echo "Processing $sample_name..."

  input_bam="${indir}/${sample_name}_nr_sorted.bam"
  filtered_bam="${outdir}/${sample_name}_sorted_q30fil_nr_sorted.bam"

  if [ ! -f "$input_bam" ]; then
    echo "Missing input BAM: $input_bam"
    continue
  fi

  echo "Filtering for MAPQ ≥ 30..."
  $samtools view -@ 8 -b -q 30 "$input_bam" > "$filtered_bam"
  $samtools index "$filtered_bam"

  echo "Finished filtering $sample_name"
done < "$samples"