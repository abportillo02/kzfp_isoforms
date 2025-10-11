#!/bin/bash
#SBATCH --job-name=merge_transcriptome
#SBATCH --output=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/merge_novel_transcriptome_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -e

# Define input files
MERGELIST="/home/abportillo/github_repo/kzfp_isoforms/scripts/mergelist_novel.txt"
REFERENCE_GTF="/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf"
REFERENCE_FASTA="/home/abportillo/github_repo/RNA_seq_Bcell/output/raw_fastq_bcell/rnaPreprocess/hg38_p14/hg38_p14.fa"

# Define output directory
OUTDIR="/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess"
MERGED_GTF="${OUTDIR}/novel_merged2.gtf"
TRANSCRIPT_FASTA="${OUTDIR}/novel_isoform_nt2.fa"
GFFCOMPARE_OUT="${OUTDIR}/gffcomp"

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

# Step 1: Merge GTFs without reference to retain all novel transcripts
echo "Merging GTFs using StringTie (without reference)..."
"$STRINGTIE" --merge -o "$MERGED_GTF" "$MERGELIST"

# Step 2: Annotate merged GTF using gffcompare
echo "Annotating merged GTF with gffcompare..."
"$GFFCOMPARE" -r "$REFERENCE_GTF" -o "$GFFCOMPARE_OUT" "$MERGED_GTF"

# Step 3: Generate transcript FASTA from annotated GTF
echo "Generating transcript FASTA using gffread..."
"$GFFREAD" "${GFFCOMPARE_OUT}.annotated.gtf" -g "$REFERENCE_FASTA" -w "$TRANSCRIPT_FASTA"

echo "Done. Output files:"
echo "  - Merged GTF: $MERGED_GTF"
echo "  - Annotated GTF: ${GFFCOMPARE_OUT}.annotated.gtf"
echo "  - Transcript FASTA: $TRANSCRIPT_FASTA"