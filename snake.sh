#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=240
#SBATCH --mem=1000

module load r r/4.3.1
module load python/3.5.1
module load fastqc/0.12.1
module load salmon/1.10.2

# Add this line when data upload is ready:
# bash scripts/setup/make_config.sh

snakemake -j 4 --latency-wait 60 --cluster "sbatch --mem=10000 -N 1 -n 12 -o ./.slurmlogs/%j.out"

snakemake --rulegraph | dot -Tsvg > dag.svg