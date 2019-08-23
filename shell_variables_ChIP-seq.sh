#!/bin/bash

# CPU threads
CPU=4
# output directory for QC
fastQC_dir="./fastQC_ChIP-seq/"
# FASTQ file directory
fastq_dir="./FASTQ_ChIP-seq/"
# file containing the list of all FASTQ files to be processed
fastq_file_list="FASTQ_files_ChIP-seq"
# FASTQ file extension
fastq_file_ext=".txt.gz$"
# file containing reference genome sequence
genome='./genome/Bacillus_subtilis_subsp._subtilis_str._168.fasta'
# BAM file directory
bam_dir="./BAM_ChIP-seq/"
# minimum MAPQ threshold for inclusion in the final BAM file
Q=10

