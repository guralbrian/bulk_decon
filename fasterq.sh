#!/bin/bash
#
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=fasterq
#SBATCH --time=2:00:00
#SBATCH --mem=5g
#SBATCH --output=./.slurmlogs/%j.out

module load sratoolkit/3.0.7

bash scripts/setup/fasterq_dump.sh