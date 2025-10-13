#!/bin/sh

#$ -m ab
#$ -m ae
#$ -j y
#$ -o out.txt
#$ -N bzip2
#$ -V


# bzip2 -dk 14893_51L-LEP-RNA_GTCCGC_L999_R1_001.fastq.bz2 

# bzip2 -dk 14895_112R-LEP-RNA_GTGGCC_L999_R1_001.fastq.bz2 

bzip2 -dk 14897_124-LEP-RNA_CGTACG_L999_R1_001.fastq.bz2 

# bzip2 -dk 14899_240L-LEP-RNA_ACTGAT_L999_R1_001.fastq.bz2

# bzip2 -dk 14901_237-LEP-RNA_CGATGT_L999_R1_001.fastq.bz2

# bzip2 -dk 14905_122L-LEP-RNA_GATCAG_L999_R1_001.fastq.bz2

# bzip2 -dk 14909_172L-LEP-RNA_AGTTCC_L999_R1_001.fastq.bz2

# bzip2 -dk 14911_191L-LEP-RNA_CCGTCC_L999_R1_001.fastq.bz2
