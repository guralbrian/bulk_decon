# 1_1_jensen_qc

# this script is meant to improve the pre-processing of the individual single cell datasets used in the pipeline
# This specifc part is meant to load each of the seurats loading in 1_jensen_load.R then filter/scale/and normalize each
# This will also implement two new methods: DoubletFinder and SoupX to remove doublets and ambient RNA contamination
# each of these will be run seperatel and saved in their own directory together 


# Load data


# load libraries
libs <- c("Seurat", "ggplot2", "DESeq2", "patchwork","SeuratDisk", "MuSiC", "reshape2",
          "tidyverse", "SingleCellExperiment","harmony", "SCpubr","shiny", 
          "AUCell", "viridis", "gplots", "scales", "ggrepel", "gridExtra", "scCustomize",
          "httr","readxl","matrixStats", "TabulaMurisSenisData") # list libraries here
lapply(libs, require, character.only = T)
source("jensen/scripts/functions/decon_all.R")

# Single cell
sn.rau <- LoadH5Seurat("jensen/data/raw/b6_sn.h5seurat")
sn.tm  <- LoadH5Seurat("jensen/data/processed/tabula_muris.h5seurat")
sn.tm$orig.ident <- paste0("tm_", sn.tm$mouse.id)
sn.tm <- subset(sn.tm, tissue_free_annotation == "Heart")

sn.wu  <- LoadH5Seurat("jensen/data/processed/single_cell/wu_2021.h5seurat")
sn.mt  <- LoadH5Seurat("jensen/data/processed/single_cell/martini_2019.h5seurat")

sn.rau$origin <- "rau"
sn.mt$origin <- "martini"


sn.rau$PercentMito <- PercentageFeatureSet(sn.rau, pattern = "^mt-")
sn.wu$PercentMito <- PercentageFeatureSet(sn.wu, pattern = "^mt-")
sn.mt$PercentMito <- PercentageFeatureSet(sn.mt, pattern = "^mt-")

# QC on each

# New ClusterSeurat includes SoupX and DoubletFinder
# and considers distribution of counts and features of each origin
sub.rau <- unique(sn.rau$orig.ident)

sn.rau <- lappy(uniqueClusterByMeta(sn.rau, sub.rau) 

# Save in new fotlder
table(sn.rau$DF.classifications_0.25_0.07_1747)
