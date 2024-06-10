libs <- c("tidyverse", "stringr", "Seurat", "SeuratDisk", "harmony", "SCpubr", 
          "ggridges", "pals", "viridis", "patchwork", "RColorBrewer") # list libraries here

libs <- c("Seurat", "SeuratDisk", "harmony", "tidyverse", "viridis", "patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

## Merge and label clusters of Rau data ##
# Get files
raw.path  <- "data/processed/single_cell/no_doublets/"
rau.files <- list.files(raw.path)[
              stringr::str_detect(list.files(raw.path), 
                             pattern = "no_doublets") 
                                  ] # Rerun last script, but save them to general folder
# Add tags for doublets and use str_detect to cut out the rest when loading
# Load and merge
rau.sn <- lapply(rau.files, function(x){LoadH5Seurat(paste0(raw.path, x))}) 
rau.sn <- merge(rau.sn[[1]], rau.sn[[2]])
rm(rau.files)

## QC and cluster ##

# Filter

rau.sn <- rau.sn |> 
            subset(nFeature_RNA   > 200     & 
                   nFeature_RNA   < 2500    &
                   nCount_RNA     > 800     &
                   PercentMito   <= 5       & 
                   DoubletScore   < 3)

# Normalize and dimensional reduction
rau.sn <- rau.sn |>
  NormalizeData(verbose = F) |>
  FindVariableFeatures(verbose = F, nfeatures = 2000) |>
  ScaleData(verbose = F) |>
  RunPCA(verbose = F) 

# Integrate replicates
rau.sn <- harmony::RunHarmony(rau.sn, 
                     group.by.vars = "orig.ident",
                     verbose = F,
                     project.dim = F)

# Find optimal PCs

# Find elbow
# Determine percent of variation associated with each PC
pct <- rau.sn[["pca"]]@stdev / sum(rau.sn[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

# Run UMAP

rau.sn <- rau.sn |>
  FindNeighbors(dims = 1:pcs, reduction = "harmony", verbose = F) |>
  FindClusters(resolution = 0.2, verbose = F) |>
  RunUMAP(dims = 1:pcs, reduction = "harmony", verbose = F)


# Save file
SaveH5Seurat(rau.sn, "data/processed/single_cell/merged_no_doublets")


# Change idents to show # of nuclei per cluster for plotting

## Save figure of UMAP with nuclei counts of each cluster in labels
temp.labels <- paste0(levels(rau.sn), " (", table(Idents(rau.sn)), " nuclei)")
names(temp.labels) <- levels(rau.sn)
rau.sn <- RenameIdents(rau.sn, temp.labels)
rau.sn$cell.type <- Idents(rau.sn)

# Color scheme
cols <-  brewer.pal(length(unique(temp.labels)), "Set3")
names(cols) <- unique(temp.labels)

if(!dir.exists("results/supp_figs/")){
  dir.create("results/supp_figs/")
}

# Save 
png(file = "results/supp_figs/cell_clusters_preannotation.png",
    width = 10, 
    height = 8,
    units = "in",
    res = 400)

# Save plot for supplement
DimPlot(rau.sn, group.by = "cell.type", cols = cols, pt.size = 1) |>
  LabelClusters(id = "cell.type", size = 6, repel = T, box = T, force = 100) +
  ggtitle(NULL) +
  theme_cowplot(font_size = 18) +
  NoLegend()

dev.off()

