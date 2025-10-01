#!/bin/bash

# Exit on error
set -e

# Define input files
MERGELIST="/home/abportillo/github_repo/kzfp_isoforms/scripts/mergelist.txt"  # List of GTF files, one per line
REFERENCE_GTF="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf"  # Reference GTF used for guided assembly
REFERENCE_FASTA="/home/abportillo/github_repo/RNA_seq_Bcell/output/raw_fastq_bcell/rnaPreprocess/hg38_p14/hg38_p14.fa"  # Reference genome FASTA

# Define output files
MERGED_GTF="merged.gtf"
TRANSCRIPT_FASTA="isoform_nt.fa"

# Define tools
STRINGTIE=/home/abportillo/.conda/envs/mamba_abner_BC/bin/stringtie
GFFREAD=/home/abportillo/.conda/envs/mamba_abner_BC/bin/gffread

# Check if required tools are available
if [ -z "$STRINGTIE" ]; then
  echo "Error: stringtie not found in PATH."
  exit 1
fi

if [ -z "$GFFREAD" ]; then
  echo "Error: gffread not found in PATH."
  exit 1
fi

# Step 1: Merge GTFs
echo "Merging GTFs using StringTie..."
$STRINGTIE --merge -G "$REFERENCE_GTF" -o "$MERGED_GTF" "$MERGELIST"

# Step 2: Generate transcript FASTA using gffread
echo "Generating transcript FASTA using gffread..."
$GFFREAD "$MERGED_GTF" -g "$REFERENCE_FASTA" -w "$TRANSCRIPT_FASTA"

echo "Done. Output files:"
echo "  - Merged GTF: $MERGED_GTF"
echo "  - Transcript FASTA: $TRANSCRIPT_FASTA"
