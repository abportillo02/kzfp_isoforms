#!/bin/sh

## get sample names and save to a txt file
ls /home/abportillo/github_repo/kzfp_isoforms/fastq_files/*.fastq | \
awk -F'_' '{print $1}' | \
sort | \
uniq > /home/abportillo/github_repo/kzfp_isoforms/scripts/sample.txt


