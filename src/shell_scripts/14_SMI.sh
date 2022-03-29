#!/usr/bin/env bash

#BSUB -J SMI
#BSUB -o logs/SMI2_%J.out
#BSUB -e logs/SMI2_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 6

module load samtools
module load R

full_path=/beevol/home/wellskri/Analysis/Lori_Sussel/Nicole_Moss/eclip_normalization

#perl peak_input_normalization_wrapper.pl cheat_manifest.txt cheat_dir
#perl peak_input_normalization_wrapper.pl manfest_file_testr2.txt testing_r2
#perl peak_input_normalization_wrapper.pl manfest_file_test_bam.txt test_bams
# perl \
#     $full_path/src/scripts/peak_input_normalization_wrapper.pl \
#     $full_path/manifest_files/manfest_file_r2.txt \
#     $full_path/data/samtools_sort_r2

perl \
    $full_path/src/scripts/peak_input_normalization_wrapper.pl \
    $full_path/manifest_files/manfest_file_test_compression_r2_perl.txt \
    $full_path/data/test_compression_r2