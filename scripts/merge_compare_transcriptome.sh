#!/bin/bash
#SBATCH --job-name=merge_compare_transcriptomes
#SBATCH --output=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/merge_compare_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -e

# Define input files
OLD_LIST="/home/abportillo/github_repo/kzfp_isoforms/scripts/old_gtf_list.txt"
YOUNG_LIST="/home/abportillo/github_repo/kzfp_isoforms/scripts/young_gtf_list.txt"
REFERENCE_GTF="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf"
REFERENCE_FASTA="/home/abportillo/github_repo/RNA_seq_Bcell/output/raw_fastq_bcell/rnaPreprocess/hg38_p14/hg38_p14.fa"

# Define output directory
OUTDIR="/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess"
MERGED_OLD="${OUTDIR}/old_merged.gtf"
MERGED_YOUNG="${OUTDIR}/young_merged.gtf"
GFFCOMPARE_OUT="${OUTDIR}/old_vs_young"
TRANSCRIPT_FASTA="${OUTDIR}/old_vs_young_transcripts.fa"

# Define tools
STRINGTIE=/home/abportillo/.conda/envs/mamba_abner_BC/bin/stringtie
GFFREAD=/home/abportillo/.conda/envs/mamba_abner_BC/bin/gffread
GFFCOMPARE=/home/abportillo/.conda/envs/mamba_abner_BC/bin/gffcompare

# Check tools
for TOOL in "$STRINGTIE" "$GFFREAD" "$GFFCOMPARE"; do
  if [ ! -x "$TOOL" ]; then
    echo "Error: $TOOL not found or not executable."
    exit 1
  fi
done

# Step 1: Merge old group
echo "Merging old group GTFs..."
"$STRINGTIE" --merge -G "$REFERENCE_GTF" -o "$MERGED_OLD" "$OLD_LIST"

# Step 2: Merge young group
echo "Merging young group GTFs..."
"$STRINGTIE" --merge -G "$REFERENCE_GTF" -o "$MERGED_YOUNG" "$YOUNG_LIST"

# Step 3: Compare merged transcriptomes
echo "Comparing old vs young merged GTFs with gffcompare..."
"$GFFCOMPARE" -r "$MERGED_YOUNG" -o "$GFFCOMPARE_OUT" "$MERGED_OLD"

# Step 4: Generate transcript FASTA from annotated comparison
echo "Generating transcript FASTA using gffread..."
"$GFFREAD" "${GFFCOMPARE_OUT}.annotated.gtf" -g "$REFERENCE_FASTA" -w "$TRANSCRIPT_FASTA"

echo "Done. Output files:"
echo "  - Old merged GTF: $MERGED_OLD"
echo "  - Young merged GTF: $MERGED_YOUNG"
echo "  - Annotated comparison GTF: ${GFFCOMPARE_OUT}.annotated.gtf"
echo "  - Transcript FASTA: $TRANSCRIPT_FASTA"