# Cardiac Cell Composition during Heart Failure

Written by [Brian Gural](https://www.linkedin.com/in/brian-gural-09bb60128/) \
README last updated on Aug 21st 2024

## Project Summary

This project is meant to measure the degree of cardiac cellular remodeling that occurs in C57BL/6 mice during heart failure, and how those compositional changes mediate bulk gene expression changes. Additionally, through a cardiomyocyte-specific knock-out model, we define the cardioprotective effect of the α1-adrenergic receptors. 

Prior research has established the role of α1-AR  [in cardiac mitocondrial function](https://pubmed.ncbi.nlm.nih.gov/35170492/) and [postinfarct remodeling](https://www.sciencedirect.com/science/article/pii/S2452302X23004187)

All data collection was performed in the labs of Drs. Brian Jensen, Christoph Rau, and Mikayla Patterson.

## Example Results

![Plot 1](https://github.com/guralbrian/bulk_decon/blob/main/results/10_plot_de/zbtb16_pik3r1.png?raw=true)

![Plot 2](https://github.com/guralbrian/bulk_decon/blob/main/results/5_findMarkers/marker_specificity.png?raw=true)

## Reproducing the analysis

This pipeline was written with SLURM-execution in mind. To run this pipeline on a SLURM-compatible HPC, simply run:

`bash snake.sh`

![DAG Plot](https://github.com/guralbrian/bulk_decon/blob/main/dag.svg)

Once the items remaining in the "To-Do" list have been finished, this should produce a copy of every figure, from the source data.

## To do:
- ~~Break scripts into functional units~~
- ~~Make pipeline executable with Snakemake~~
- Containerize with Singularity or version control with .Renv
- Download SRA files
- ~~Add .fastq pseudoalignment~~



