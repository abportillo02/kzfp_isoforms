#!/bin/bash
#SBATCH --job-name=salmon_quant_novel
#SBATCH --output=/home/abportillo/github_repo/kzfp_isoforms/salmon_quant_novel_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -e

# Define paths
TRANSCRIPT_FASTA="/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/novel_isoform_nt3.fa"
INDEX_DIR="/home/abportillo/github_repo/kzfp_isoforms/salmon_index_novel3"
READ_DIR="/home/abportillo/github_repo/kzfp_isoforms/fastq_files"
OUTPUT_DIR="/home/abportillo/github_repo/kzfp_isoforms/salmon_quant_novel3"

# Load Salmon if needed (or use full path)
SALMON=/home/abportillo/.conda/envs/mamba_abner_BC/bin/salmon

if [ -z "$SALMON" ]; then
  echo "Error: Salmon not found in PATH."
  exit 1
fi

# Create output directories
mkdir -p "$INDEX_DIR"
mkdir -p "$OUTPUT_DIR"

# Step 1: Build Salmon index
echo "Building Salmon index..."
$SALMON index -t "$TRANSCRIPT_FASTA" -i "$INDEX_DIR"

# Step 2: Quantify each sample
for fq in "$READ_DIR"/*.fastq; do
  sample=$(basename "$fq" .fastq)
  echo "Quantifying $sample..."
  $SALMON quant -i "$INDEX_DIR" -l A \
    -r "$fq" \
    -p 8 \
    --validateMappings \
    --seqBias --gcBias \
    -o "$OUTPUT_DIR/${sample}_quant"
done

echo "Salmon quantification complete."