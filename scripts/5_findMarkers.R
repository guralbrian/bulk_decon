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
Idents(sn.mark)[which(Idents(sn.mark) == "4")] <- "0"

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
    #dplyr::filter(summary.logFC >= 1 |
    #                summary.logFC <= -1) |>
    arrange(FDR) |>
    #slice_head(n = 10) |>
    dplyr::select(p.value, FDR, summary.logFC) |>
    mutate(celltype = type) |> 
    rownames_to_column(var = "gene")
  return(marker)}

# Collapse list of data frames
all.markers <- lapply(levels(Idents(sn.mark)), function(x){.getMarkers(x)}) |>
  purrr::reduce(full_join)

# get the top markers
top.markers <- all.markers |> 
  dplyr::filter(summary.logFC >= 1 |
                  summary.logFC <= -1) |>
  arrange(FDR) |>
  slice_head(n = 10) 
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
                "VSMC/Pericytes",
                "Endothelial Cells",
                "Macrophages",
                "SMC")

# Subset to the high-confidence clusters
# Rename the clusters to match the cell types
sn.mark <- sn |>
  subset(idents = seq(0,length(cell.types)-1,1)) 
names(cell.types) <- levels(sn.mark)
  sn.mark <- RenameIdents(sn.mark, cell.types)
sn.mark$cell.type <- Idents(sn.mark)

# Save UMAP plot
# Color scheme
cols <-  brewer.pal(length(cell.types), "Set2")
names(cols) <- cell.types

# Save 
png(file = "results/5_findMarkers/cell_clusters.png",
    width = 1000, 
    height = 800,
    units = "px",
    res = 100)

DimPlot(sn.mark, group.by = "cell.type", cols = cols, pt.size = 1) |>
  LabelClusters(id = "cell.type", size = 6, repel = T, box = T, force = 100) +
  ggtitle(NULL) +
  theme_cowplot(font_size = 18) +
  NoLegend()


dev.off()


# Save markers, add cell types and removed unused cell types
top.markers <- top.markers |>
  mutate(annotation = cell.types[as.character(celltype)]) |>
  subset(!is.na(annotation))

# Save markers, add cell types and removed unused cell types
all.markers <- all.markers |>
  mutate(annotation = cell.types[as.character(celltype)]) |>
  subset(!is.na(annotation))

# Save the data
write.csv(all.markers, "data/processed/single_cell/all_markers.csv", row.names = F)
write.csv(top.markers, "data/processed/single_cell/cluster_markers.csv", row.names = F)
SaveH5Seurat(sn.mark, "data/processed/single_cell/celltype_labeled",  overwrite = T)
