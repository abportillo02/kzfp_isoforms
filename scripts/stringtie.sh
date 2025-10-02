#!/bin/bash

# Check for sample list
if [ -z "$1" ]; then
  echo "Usage: $0 <sample_list.txt>"
  exit 1
fi

samples="$1"
bamdir=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/filtered_bams
outdir=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/stringtie_ballgown
# gtf=/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/teAnno_round3/filtered_gencode_v46_chr_patch_hapl_scaff_annotation.gtf
stringtie=/home/abportillo/.conda/envs/mamba_abner_BC/bin/stringtie

mkdir -p "$outdir"

while IFS= read -r sample_name; do
  echo "Creating job for $sample_name..."

  job_script="${outdir}/${sample_name}_stringtie_job.sh"
  input_bam="${bamdir}/${sample_name}_nr_sorted.bam"
  sample_outdir="${outdir}/ballgown/${sample_name}"

  mkdir -p "$sample_outdir"

  cat <<EOF>> "$job_script"
#!/bin/bash
#SBATCH --job-name=StringTie_${sample_name}
#SBATCH --output=${outdir}/${sample_name}_stringtie_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=64G
#SBATCH --time=24:00:00

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

$stringtie "$input_bam" -e -B -p 8 \
  # -G "$gtf" \
  -o "${sample_outdir}/${sample_name}_novel_ballgown.gtf" \
  -A "${sample_outdir}/counts.txt"
EOF

  sbatch "$job_script"
done < "$samples"