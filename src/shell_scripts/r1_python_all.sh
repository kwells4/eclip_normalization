#!/usr/bin/env bash

#BSUB -J SMI
#BSUB -o logs/SMI2_%J.out
#BSUB -e logs/SMI2_%J.err
#BSUB -R "select[mem>64] rusage[mem=64] " 
#BSUB -q rna
#BSUB -n 1

module load samtools
module load R

full_path=/beevol/home/wellskri/Analysis/Lori_Sussel/Nicole_Moss/eclip_normalization

# python \
#     $full_path/src/scripts/peak_input_normalization.py \
#     -m $full_path/manifest_files/manfest_file_r1_python.txt \
#     -o $full_path/python_results_r1_perl_logic \
#     -f perl_script

python \
    $full_path/src/scripts/peak_input_normalization.py \
    -m $full_path/manifest_files/manfest_file_r1_python.txt \
    -o $full_path/python_results_r1 \
    -f default