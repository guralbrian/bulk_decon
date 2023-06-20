# 1_1_jensen_qc

# this script is meant to improve the pre-processing of the individual single cell datasets used in the pipeline
# This specifc part is meant to load each of the seurats loading in 1_jensen_load.R then filter/scale/and normalize each
# This will also implement two new methods: DoubletFinder and SoupX to remove doublets and ambient RNA contamination
# each of these will be run seperatel and saved in their own directory together 

# load libraries
libs <- c("Seurat", "ggplot2", "DESeq2", "patchwork","SeuratDisk", "reshape2",
          "tidyverse", "SingleCellExperiment", "harmony", "SCpubr", 
          "AUCell", "viridis", "gplots", "scales", "ggrepel", "gridExtra", "scCustomize",
          "httr","matrixStats", "scran", "scuttle", "scater", "DropletUtils", "scDblFinder") # list libraries here

lapply(libs, require, character.only = T)
source("jensen/scripts/functions/decon_all.R")

sn.datasets <- list.files("jensen/data/raw/single_cell/")

DropletsAllDatasets <- function(sn.dataset){

# name of dataset
dataset <- sn.dataset

# path to dataset
file <- paste0("jensen/data/raw/single_cell/", dataset)

# make a directory to store the results
dir.path <- paste0("jensen/results/doublets_ambient_rna/", strsplit(dataset, "[.]")[[1]][[1]])
if(!dir.exists(dir.path)) {
  dir.create(dir.path)
}
# load the data
sn <- LoadH5Seurat(file)

# QC on each sample

# New ClusterSeurat includes SoupX and DoubletFinder
# and considers distribution of counts and features of each origin
sub.list <- unique(sn$orig.ident)
result_list <-lapply(sub.list, RemoveDoublets, seurat = sn, directory = dir.path)
result <- Reduce(
  f = function(x, y) {merge(x, y, merge.data = FALSE)},
  x = result_list # list of Seurat objects
)
print(paste0("reduced seurat list for ", dataset))
SaveH5Seurat(result, paste0("jensen/data/processed/single_cell/no_doublets/", dataset))
}

lapply(sn.datasets, function(x){DropletsAllDatasets(x)})
### you left off here
# summary of last things #

# looked at https://doi.org/10.1186/s13059-023-02978-x
# looked for other droplet detection methods that work with current data (counts matrix)
# following this tutorial https://rockefelleruniversity.github.io/scRNA-seq/exercises/answers/exercise2_answers.html

# need to do wu dataset still
