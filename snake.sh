#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=240
#SBATCH --mem=1000

module load r r/4.3.1
module load python/3.9.6
module load fastqc/0.12.1
module load salmon/1.10.2
module load sratoolkit/3.0.7

# Add this line when data upload is ready:
# bash scripts/setup/make_config.sh

snakemake -j 36 --latency-wait 60 --cluster "sbatch --mem=128000 -N 1 -n 12 -o ./.slurmlogs/%j.out" 

snakemake --rulegraph | dot -Tsvg > dag.svg