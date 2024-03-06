#!/bin/bash

# Save the current directory
parent_dir=$(pwd)

# Navigate to the Directory and List Contents
cd data/raw/fastq/

# List Directory Names
dirs=$(ls -d */) 

# Replace Hyphens with Underscores and Remove Trailing Slashes
samples=$(echo $dirs | tr '-' '_' | tr -d '/')

# Format as JSON Array
json_samples=$(echo $samples | sed 's/ /", "/g')

# Create the JSON Structure
json_content="{\"samples_fastq\": [\"$json_samples\"]}"

# Navigate back to the parent directory
cd $parent_dir

# Write to fastq_config.json in the parent directory
echo $json_content > scripts/setup/fastq_config.json

