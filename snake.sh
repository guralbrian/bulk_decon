#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=240
#SBATCH --mem=1000

module load r r/4.2.1
module load python
module load fastqc
module load fastqc
module load salmon

# Add this line when data upload is ready:
# bash scripts/setup/make_config.sh

snakemake -j 4 --latency-wait 30 --cluster "sbatch --mem=10000 -N 1 -n 12 -o ./.slurmlogs/%j.out"

snakemake --dag | dot -Tsvg > dag.svg