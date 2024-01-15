# Find markers + annotate clusters

libs <- c("tidyverse", "Seurat", "SeuratDisk",  "Biobase", 
          "reshape2", "SingleCellExperiment") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)


sn <- LoadH5Seurat("data/processed/single_cell/merged_no_doublets.h5seurat")

# Load whole bulk RNAseq
bulk.all <- read.csv("data/processed/bulk/all_counts.csv")


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
    dplyr::filter(summary.logFC >= 1 |
                    summary.logFC <= -1) |>
    arrange(FDR) |>
    slice_head(n = 10) |>
    dplyr::select(p.value, FDR, summary.logFC) |>
    mutate(celltype = type) |> 
    rownames_to_column(var = "gene")
  return(marker)}

# Collapse list of data frames
all.markers <- lapply(levels(Idents(sn.mark)), function(x){.getMarkers(x)}) |>
  purrr::reduce(full_join)


gc()

# Genes here are printed and manually entered into ToppGene
# https://toppgene.cchmc.org/
genes <- function(x){ all.markers |> 
    subset(celltype == x) |> 
    pull(gene) |> 
    as.data.frame() |> 
    print(row.names = F, print.keys = F)
}

cell.types <- c("Endothelial Cells",
                "Cardiomyocytes",
                "Fibroblast",
                "VSMC/Pericytes",
                "Endothelial Cells 2",
                "Macrophages",
                "SMC")
#"Cardiac Neuron", # Neuron is the last cell type to be identified with confidence
#"?Matrix FB, Lung?")
#"unclear",#"??EC, CA??"
#"unclear")#"???")

# Subset to the high-confidence clusters
# Rename the clusters to match the cell types
sn.mark <- sn |>
  subset(idents = seq(0,length(cell.types)-1,1)) 
#subset(idents = c(0,1,2,4,5))
names(cell.types) <- levels(sn.mark)
sn.mark <- RenameIdents(sn.mark, cell.types)

row.names(bulk.es.exp) <- all.markers$gene


# Save markers, add cell types and removed unused cell types
all.markers <- all.markers |>
  mutate(annotation = cell.types[as.character(celltype)]) |>
  subset(!is.na(annotation))

# Save the data
write.csv(all.markers, "data/processed/single_cell/cluster_markers.csv", row.names = F)
SaveH5Seurat(sn.mark, "data/processed/single_cell/celltype_labeled")
