#!/bin/bash
#SBATCH --time=725:00:00
#SBATCH --job-name="A long job"
#SBATCH --output=long-job.out

Rscript sim_main_hd.R > output_funksim_$1_$2_$3_$4_$5_$6.txt p=$1 s=$2 rho=$3 design=$4 iter=$5 method=$6
