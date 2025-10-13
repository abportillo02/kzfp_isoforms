#!/bin/bash
#SBATCH --job-name=generate_fasta
#SBATCH --output=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/generate_fasta_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -e

# Define input files
# MERGELIST="/home/abportillo/github_repo/kzfp_isoforms/scripts/mergelist.txt"
# REFERENCE_GTF="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf"
REFERENCE_FASTA="/home/abportillo/github_repo/RNA_seq_Bcell/output/raw_fastq_bcell/rnaPreprocess/hg38_p14/hg38_p14.fa"

# Define output directory
OUTDIR="/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess"
MERGED_GTF="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/te_function/rna_seq_output/young_old_LEP_R1_gemini/stringtie_merge/old_vs_young.annotated.gtf"
TRANSCRIPT_FASTA="${OUTDIR}/qianhui_output.fa"

# Define tools
STRINGTIE=/home/abportillo/.conda/envs/mamba_abner_BC/bin/stringtie
GFFREAD=/home/abportillo/.conda/envs/mamba_abner_BC/bin/gffread

# Check tools
if [ ! -x "$STRINGTIE" ]; then
  echo "Error: stringtie not found or not executable."
  exit 1
fi

if [ ! -x "$GFFREAD" ]; then
  echo "Error: gffread not found or not executable."
  exit 1
fi

# # Step 1: Merge GTFs
# echo "Merging GTFs using StringTie..."
# "$STRINGTIE" --merge -G "$REFERENCE_GTF" -o "$MERGED_GTF" "$MERGELIST"

#Step 2: Generate transcript FASTA
echo "Generating transcript FASTA using gffread..."
"$GFFREAD" "$MERGED_GTF" -g "$REFERENCE_FASTA" -w "$TRANSCRIPT_FASTA"

echo "Done. Output files:"
echo "  - Merged GTF: $MERGED_GTF"
echo "  - Transcript FASTA: $TRANSCRIPT_FASTA"