# Find markers + annotate clusters

libs <- c("tidyverse", "Seurat", "SeuratDisk",  "Biobase", 
          "reshape2", "SingleCellExperiment", "RColorBrewer", "cowplot") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)


sn <- LoadH5Seurat("data/processed/single_cell/merged_no_doublets.h5seurat")

# Load whole bulk RNAseq
bulk.all <- read.csv("data/processed/bulk/all_counts.csv", row.names = 1)


# Find markers for sn clusters, annotate to cell types
sn.mark <- sn |>
  subset(features = rownames(sn)[rownames(sn) %in% rownames(bulk.all)])

sn.sce <- sn.mark |> 
  as.SingleCellExperiment(assay = "RNA") |>
  as("SummarizedExperiment")

# Add gene names
rownames(sn.sce) <- rownames(sn.mark@assays$RNA@counts)

# Get markers based on expression ratios between clusters and 1vAll comparisons

scn.markers  <- scran::findMarkers(sn.sce, groups = Idents(sn.mark), pval.type = "some", assay.type = "counts")

.getMarkers <- function(type){
  marker <- scn.markers@listData[[type]] |>
    as.data.frame() %>%
    dplyr::select(p.value, FDR, summary.logFC) |>
    mutate(celltype = type) |> 
    rownames_to_column(var = "gene")
  return(marker)}

# Collapse list of data frames
all.markers <- lapply(levels(Idents(sn.mark)), function(x){.getMarkers(x)}) |>
  purrr::reduce(full_join)

# get the top markers
top.markers <- all.markers |> 
  group_by(celltype) |> 
  dplyr::filter(summary.logFC >= 1) |>#|
                 # summary.logFC <= -1) |>
  arrange(p.value) |>
  slice_head(n = 15) 
gc()

# Genes here are printed and manually entered into ToppGene
# https://toppgene.cchmc.org/
genes <- function(x){ top.markers |> 
    subset(celltype == x) |> 
    pull(gene) |> 
    as.data.frame() |> 
    print(row.names = F, print.keys = F)
}

cell.types <- c("Endothelial Cells",
                "Cardiomyocytes",
                "Fibroblast",
                "Pericytes/SMC",
                "Macrophage",
                "Pericytes/SMC")
                #"Endothelial Cells","Cardiac Neuron","Macrophage"

# Subset to the high-confidence clusters and redo marker ID
# Rename the clusters to match the cell types
sn.mark <- sn |>
  subset(idents = seq(0,length(cell.types)-1,1)) 
names(cell.types) <- levels(sn.mark)
sn.mark <- RenameIdents(sn.mark, cell.types)
sn.mark$cell.type <- Idents(sn.mark)


sn.sce <- sn.mark |> 
  as.SingleCellExperiment(assay = "RNA") |>
  as("SummarizedExperiment")

# Add gene names
rownames(sn.sce) <- rownames(sn.mark@assays$RNA@counts)

# Rerun marker selection with all pval type
scn.markers  <- scran::findMarkers(sn.sce, groups = Idents(sn.mark), pval.type = "all", assay.type = "counts")

# Collapse list of data frames
all.markers <- lapply(levels(Idents(sn.mark)), function(x){.getMarkers(x)}) |>
  purrr::reduce(full_join)

# get the top markers
top.markers <- all.markers |> 
  group_by(celltype) |> 
  dplyr::filter(summary.logFC >= 1
              | summary.logFC <= 1) |>
  arrange(p.value) |>
  slice_head(n = 15) 

# Save UMAP plot
# Color scheme
cols <-  brewer.pal(length(unique(cell.types)), "Set2")
names(cols) <- unique(cell.types)

# Save 
png(file = "results/5_findMarkers/cell_clusters.png",
    width = 10, 
    height = 8,
    units = "in",
    res = 400)

DimPlot(sn.mark, group.by = "cell.type", cols = cols, pt.size = 1) |>
  LabelClusters(id = "cell.type", size = 6, repel = T, box = T, force = 100) +
  ggtitle(NULL) +
  theme_cowplot(font_size = 18) +
  NoLegend()


dev.off()

# Save the data
write.csv(all.markers, "data/processed/single_cell/all_markers.csv", row.names = F)
write.csv(top.markers, "data/processed/single_cell/cluster_markers.csv", row.names = F)
SaveH5Seurat(sn.mark, "data/processed/single_cell/celltype_labeled",  overwrite = T)
