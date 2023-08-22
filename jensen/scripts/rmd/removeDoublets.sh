#!/bin/bash


#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 20:00:00
#SBATCH --mem=32g
#SBATCH -n 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bgural@email.unc.edu

module load r/4.2.1

Rscript by_date/08142023_doublet_ambient.R
