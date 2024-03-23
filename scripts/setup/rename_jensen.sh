#!/bin/bash

# Navigate to the directory containing the fastq files
cd data/raw/fastq/

# Loop over files that do not start with 'B6_' or 'multiqc'
for file in *; do
    if [[ ! $file =~ ^B6_ ]] && [[ ! $file =~ ^multiqc ]]; then
        # Extract the prefix up to the last underscore
        prefix=$(echo $file | rev | cut -d'_' -f2- | rev)
        
        # Rename the file to end with '.fastq.gz' instead of '.fq.gz'
        new_name=$(echo $file | sed 's/\.fq\.gz$/.fastq.gz/')

        # Create a directory for this prefix if it does not already exist
        mkdir -p "$prefix"

        # Move and rename the file into its corresponding directory
        mv "$file" "${prefix}/${new_name}"
    fi
done