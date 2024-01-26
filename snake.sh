#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=1440
#SBATCH --mem=1000

module load r r/4.2.1
module load python

snakemake -j 2 --latency-wait 30 --cluster "sbatch --mem=16000 -N 1 -n 1 -o ./.slurmlogs/%j.out"

snakemake --dag | dot -Tsvg > dag.svg