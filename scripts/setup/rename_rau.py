#!/usr/bin/env python3
# This script is just meant to rename the raw .fastq files supplied to me by christoph

import os
import re

def rename_files(directory_path):
    # Pattern to match strings starting with 'B6' and ending with 'RNA_{any number}'
    pattern = re.compile(r"(B6.*RNA_\d+)_(S\d+_L002_R)([12])_001\.fastq\.gz")

    for root, dirs, files in os.walk(directory_path):
        for filename in files:
            if filename.endswith(".fastq.gz"):
                match = pattern.match(filename)
                if match:
                    new_filename = "{}_{}.fastq.gz".format(match.group(1), match.group(3))
                    old_file_path = os.path.join(root, filename)
                    new_file_path = os.path.join(root, new_filename)
                    os.rename(old_file_path, new_file_path)
                    print("Renamed '{}' to '{}'".format(filename, new_filename))

# Base directory containing all subdirectories with .fastq.gz files
base_directory_path = "data/raw/fastq"
rename_files(base_directory_path)
