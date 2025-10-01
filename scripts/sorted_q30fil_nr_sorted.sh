#!/bin/bash

# Check if input file with sample names is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <sample_list.txt>"
  exit 1
fi

samples="$1"
outdir=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/Star_files
samtools=/home/abportillo/.conda/envs/mamba_abner_BC/bin/samtools

while IFS= read -r sample_name; do
  echo "Processing $sample_name..."

  input_bam="${outdir}/${sample_name}_Aligned.out.bam"
  sorted_bam="${outdir}/${sample_name}_sorted.bam"
  filtered_bam="${outdir}/${sample_name}_q30.bam"
  final_bam="${outdir}/${sample_name}_q30_sorted.bam"

  # Step 1: Sort original BAM
  if [ ! -f "$input_bam" ]; then
    echo "Missing input BAM: $input_bam"
    continue
  fi

  echo "Sorting original BAM..."
  $samtools sort -@ 8 -O bam -o "$sorted_bam" "$input_bam"
  $samtools index "$sorted_bam"

  # Step 2: Filter for MAPQ ≥ 30
  echo "Filtering for MAPQ ≥ 30..."
  $samtools view -@ 8 -b -q 30 "$sorted_bam" > "$filtered_bam"

  # Step 3: Sort filtered BAM
  echo "Sorting filtered BAM..."
  $samtools sort -@ 8 -O bam -o "$final_bam" "$filtered_bam"
  $samtools index "$final_bam"

  echo "Finished processing $sample_name"
done < "$samples"

  SORT_ORDER=coordinate