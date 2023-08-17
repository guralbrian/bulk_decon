# Load libraries

libs <- c("Seurat", "ggplot2", "DESeq2", "patchwork","SeuratDisk", "MuSiC", "reshape2",
          "tidyverse", "SingleCellExperiment", "harmony", "SCpubr", 
          "AUCell", "viridis", "gplots", "scales", "ggrepel", "gridExtra", "scCustomize",
          "httr","matrixStats", "scran", "scuttle", "scater", "DropletUtils", "scDblFinder") # list libraries here

lapply(libs, require, character.only = T)
source("jensen/scripts/functions/decon_all.R")

# load individual seurats
doublet.path <- "jensen/data/processed/single_cell/no_doublets/06202023/"
datasets <- list.files(doublet.path)

sn.list <- lapply(datasets, function(x){LoadH5Seurat(paste0(doublet.path, x))})

# filter each of them by quantile values
sn.list <- lapply(sn.list, function(x){FilterByQuantile(x, pt.remove = 0.1)})

# Join into one
sn.all <- merge(sn.list[[1]], sn.list[[2]]) |>
            merge(sn.list[[3]]) |>
            merge(sn.list[[4]])

# Annotate to cell types

markers.broad <- read.csv("jensen/data/processed/external/mclellan_2020/mclellan_cell_markers_broad.csv")

# Split 'markers' into separate data frames for each unique subcluster
subclusters <- split(markers.broad, markers.broad$cluster)
# Extract the gene names from each subcluster data frame
subcluster_genes <- lapply(subclusters, function(x) x$gene)
# Create a list with named elements corresponding to each subcluster
geneSets <- setNames(subcluster_genes, names(subclusters))
# make expression matrix of single cell
my.expr <-  sn.all@assays$RNA@counts
# Find  AUC for each cell by each type
cells_AUC <- AUCell_run(my.expr, geneSets, )
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE)

# get the thresholds and report when/how they're passdd
selectedThresholds <- getThresholdSelected(cells_assignment)

auc.val <- cells_AUC@assays@data@listData$AUC |>
  as.data.frame() |>
  t()

# Apply process_cell_id function to each row of auc.val using lapply()
broad_cell_types_list <- lapply(1:nrow(auc.val), function(i) processCellId(auc.val[i, ]))

# Convert the list to a data frame
broad.cell.types <- data.frame(BroadCellType = unlist(broad_cell_types_list), stringsAsFactors = FALSE)
rownames(broad.cell.types) <- rownames(auc.val)

sn.all <- AddMetaData(sn.all, broad.cell.types)

SaveH5Seurat(sn.all, "jensen/data/processed/single_cell/merged_datasets/no_doublets_2")

