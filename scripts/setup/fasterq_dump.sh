#!/bin/bash

# This script is used to generate fastq files from the raw SRA data
# It iteratively loops through each SRA, converting with fasterq-dump

# Define the path to the SRA files
SRA_PATH="./data/raw/sra_files"

# Define the path to the output fastq files
FASTQ_PATH="./data/raw/fastq_sra"

# Make vector of directory names in SRA_PATH
SRA_DIRS=$(ls $SRA_PATH)

# Loop through each directory
for dir in $SRA_DIRS
do
    # Define the full path to the SRA directory
    SRA_DIR=$SRA_PATH/$dir

    # Define the full path to the output fastq directory
    FASTQ_DIR=$FASTQ_PATH/$dir

    # Check if the fastq directory already exists
    if [ -f $FASTQ_DIR ]
    then
        echo "Fastq output directory already exists: $FASTQ_DIR"
    else
        # Convert the SRA file to fastq
        fasterq-dump $SRA_DIR -O $FASTQ_DIR
    fi
done