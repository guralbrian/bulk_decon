#!/usr/bin/env python3

import os
import re
import shutil
import zipfile

# Define the path to the zip file and the extraction directory
zip_path = 'data/raw/fastq/RNAseq.zip'
extract_dir = os.path.dirname(zip_path)

# Extract the zip file
with zipfile.ZipFile(zip_path, 'r') as zip_ref:
    zip_ref.extractall(extract_dir)

# Loop through the files in the extraction directory
for filename in os.listdir(extract_dir):
    if filename.endswith('.fq.gz'):
        # Identify the unique pattern preceding '_1.fq.gz' or '_2.fq.gz'
        match = re.match(r'(.+?)_[12]\.fq\.gz', filename)
        if match:
            unique_pattern = match.group(1)
            # Create a directory for the unique pattern if it doesn't exist
            dir_path = os.path.join(extract_dir, unique_pattern)
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
            # Move the file to its corresponding directory
            shutil.move(os.path.join(extract_dir, filename), dir_path)

print("Files have been organized into directories based on their unique patterns.")
