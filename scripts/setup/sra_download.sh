#!/bin/bash

# Make and set the directory for the files to load into
output_dir="data/raw/sra_files"
mkdir -p $output_dir
cd $output_dir

# Study ID
study_id="SRP513262"

# Download the SRA files
prefetch -v $study_id

