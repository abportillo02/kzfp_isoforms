#!/bin/sh

# Check if the input file is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

samples=$1

# Check if the input file exists
if [ ! -f "$samples" ]; then
  echo "Input file not found!"
  exit 1
fi

# Read each line from the input file and create a job script for each sample
while IFS= read -r sample_name; do
  echo "Creating job script for sample: $sample_name"

  datapath_kzfp_isoform=/home/abportillo/github_repo/kzfp_isoforms/fastq_files
  mkdir -p /home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess
  mkdir -p /home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/Star_files_2
  outdir=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/rnaPreprocess/Star_files_2
  script_path=${outdir}/${sample_name}_star_job.sh

  java=/home/abportillo/.conda/envs/mamba_abner_BC/bin/java
  bamCoverage=/home/abportillo/.conda/envs/mamba_abner_BC/bin/bamCoverage
  samtools=/home/abportillo/.conda/envs/mamba_abner_BC/bin/samtools
  STAR=/home/abportillo/.conda/envs/mamba_abner_BC/bin/STAR
  picard=/home/abportillo/.conda/envs/mamba_abner_BC/bin/picard
  wigToBigWig=/home/abportillo/.conda/envs/mamba_abner_BC/bin/wigToBigWig

  cat <<EOF > "$script_path"
#!/bin/bash

#SBATCH --job-name=RNA_hg38_p14_2passStar
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abportillo@coh.org
#SBATCH -n 16
#SBATCH -N 1-4
#SBATCH -p all
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=/home/abportillo/github_repo/kzfp_isoforms/fastq_files/14897_124-LEP-RNA_CGTACG_L999_R1_001_RNA_hg38_p14_2passStar_%j.log

source /home/abportillo/.bashrc
conda activate /home/abportillo/.conda/envs/mamba_abner_BC

mkdir -p ${outdir}/fastqc_out
module load FastQC/0.11.8
fastqc -t 8 -o ${outdir}/fastqc_out ${datapath_kzfp_isoform}/${sample_name}.fastq
module unload FastQC/0.11.8

${STAR} --genomeDir /home/abportillo/github_repo/RNA_seq_Bcell/output/raw_fastq_bcell/rnaPreprocess/hg38_p14/STAR_hg38_p14_geneCodeGTF_filter \
--readFilesIn ${datapath_kzfp_isoform}/${sample_name}.fastq \
--runThreadN 8 \
--twopassMode None \
--outFilterMultimapNmax 10 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 3 \
--outSJfilterReads Unique \
--outSJfilterOverhangMin 10 10 10 10 \
--outFileNamePrefix ${outdir}/${sample_name}_ \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMstrandField intronMotif \
--outSAMunmapped Within \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
--chimOutJunctionFormat 1 \
--chimMainSegmentMultNmax 1 \
--outSAMattributes NH HI AS nM NM ch



${samtools} sort -@ 8 -O bam -o ${outdir}/${sample_name}_sorted.bam ${outdir}/${sample_name}_Aligned.out.bam
${samtools} index ${outdir}/${sample_name}_sorted.bam
rm ${outdir}/${sample_name}_Aligned.out.bam

${picard} AddOrReplaceReadGroups \
-I ${outdir}/${sample_name}_sorted.bam \
-O ${outdir}/${sample_name}_rg_sorted.bam \
--RGID ${sample_name} \
--RGLB default_lib \
--RGPL ILLUMINA \
--RGPU unknown \
--RGSM ${sample_name} \
--SORT_ORDER coordinate

${java} -Djava.io.tmpdir=/home/abportillo/github_repo/RNA_seq_Bcell/scripts/raw_fastq_bcell/rnaPreprocess/temp \
-jar /home/abportillo/.conda/envs/mamba_abner_BC/share/picard-3.3.0-0/picard.jar MarkDuplicates \
--INPUT ${outdir}/${sample_name}_rg_sorted.bam \
--OUTPUT ${outdir}/${sample_name}_nr_sorted.bam \
--REMOVE_DUPLICATES true \
--READ_NAME_REGEX null \
--METRICS_FILE ${outdir}/${sample_name}_picardStats.txt

${samtools} sort -@ 8 -O bam -o ${outdir}/${sample_name}_sorted_nr_sorted.bam ${outdir}/${sample_name}_nr_sorted.bam
${samtools} index -@ 8 ${outdir}/${sample_name}_sorted_nr_sorted.bam
EOF

  # Submit the job
  sbatch "$script_path"
done < "$samples"