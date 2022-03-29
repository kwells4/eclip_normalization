#!/usr/bin/env bash

#BSUB -J SMI
#BSUB -o logs/SMI2_%J.out
#BSUB -e logs/SMI2_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 6

module load samtools
module load R

#perl peak_input_normalization_wrapper.pl cheat_manifest.txt cheat_dir
#perl peak_input_normalization_wrapper.pl manfest_file_testr2.txt testing_r2
perl peak_input_normalization_wrapper.pl manfest_file_test_bam.txt test_bams