#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --mem=4g
#SBATCH -n 1
#SBATCH -o ./.slurmlogs/%j.out

module load r r/4.3.1

# Run the R script with the model index and chunk index as arguments
Rscript scripts/12_02_simulate_deseq2.R $1 $2