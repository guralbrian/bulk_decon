# This script is part of the initial processing pipeline for Christoph's cell fraction RNAseq data
# In this script, we will compile the quant.sf files into a single count matrix
# Then, counts will be converted from transcript- to gene-leve anotations with tximport

# Load libs
libs <- c("tidyverse", "tximport") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# List sample quant.sf files
files <- list.files("data/raw/fastq")
files <- files[str_detect(files, "B6_")]
sample.files <- file.path("data/raw/fastq", files, "quant.sf")
names(sample.files) <- files

all(file.exists(sample.files))
# Make transcript + gene dict w/ tx2gene
vignette("tximportData")

txi <- tximport(sample.files, type = "salmon")
