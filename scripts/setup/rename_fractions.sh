#!/bin/bash

# First, rename all files
find data/raw/fastq/ -type f -name "*-*" | while IFS= read -r file; do
    dir=$(dirname "$file")
    base=$(basename "$file")
    newname=$(echo "$base" | sed 's/-/_/g')
    mv "$file" "$dir/$newname"
done

# Then, rename directories
find data/raw/fastq/ -type d -depth -name "*-*" | while IFS= read -r dir; do
    newdir=$(echo "$dir" | sed 's/-/_/g')
    mv "$dir" "$newdir"
done
