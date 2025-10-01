#!/bin/sh

# Get sample names from .fastq files only (excluding .bz2)
ls /home/abportillo/github_repo/kzfp_isoforms/fastq_files/*.fastq | \
grep -v '.bz2' | \
awk -F'/' '{print $NF}' | \
sed 's/\.fastq$//' | \
sort | uniq > /home/abportillo/github_repo/kzfp_isoforms/scripts/sample.txt

