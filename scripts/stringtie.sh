#!/bin/bash

# Check if input file with sample names is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <sample_list.txt>"
  exit 1
fi

samples="$1"
bamdir=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/filtered_bams
outdir=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/stringtie_ballgown
gtf=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf
stringtie=/home/abportillo/.conda/envs/mamba_abner_BC/bin/stringtie

# Create output directory if it doesn't exist
mkdir -p "$outdir"

while IFS= read -r sample_name; do
  echo "Running StringTie for $sample_name..."

  input_bam="${bamdir}/${sample_name}_sorted_q30fil_nr_sorted.bam"
  output_gtf="${outdir}/${sample_name}_Gencode_transcripts_ballgown.gtf"
  sample_outdir="${outdir}/ballgown/${sample_name}"

  if [ ! -f "$input_bam" ]; then
    echo "Missing BAM file: $input_bam"
    continue
  fi

  mkdir -p "$sample_outdir"

  $stringtie "$input_bam" -e -B -p 8 \
    -G "$gtf" \
    -o "$output_gtf" \
    -A "${sample_outdir}/counts.txt"

  echo "Finished $sample_name"
done < "$samples"