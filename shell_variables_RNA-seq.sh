#!/bin/bash

# CPU threads
CPU=4
# output directory for QC
fastQC_dir="./fastQC_RNA-seq/"
# FASTQ file directory
fastq_dir="./FASTQ_RNA-seq/"
# file containing the list of all FASTQ files to be processed
fastq_file_list="FASTQ_files_RNA-seq"
# FASTQ file extension
fastq_file_ext=".txt.gz$"
# file containing reference genome sequence
genome='./genome/Bacillus_subtilis_subsp._subtilis_str._168.fasta'
# BAM file directory
bam_dir="./BAM_RNA-seq/"
# minimum MAPQ threshold for inclusion in the final BAM file
Q=10
# output directory for genome coverage data
coverage_dir="./coverage_RNA-seq/"
# file containing the list of all BAM files to be processed
bam_file_list="BAM_files_RNA-seq"
# string contained in KO sample file names
ko="KO"
# string contained in WT sample file names
wt="WT"
