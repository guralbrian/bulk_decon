libs <- c("tidyverse", "stringr", "Seurat", "SeuratDisk", "harmony", "SCpubr", 
          "ggridges", "pals", "viridis", "patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

## Merge and label clusters of Rau data ##
# Get files
raw.path  <- "data/processed/single_cell/"
rau.files <- list.files(raw.path)[
              stringr::str_detect(list.files(raw.path), 
                             pattern = "rau") 
                                  ]
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
SaveH5Seurat(rau.sn, "data/processed/single_cell/merges/rau_patterson/09132023/doublets_below_3")

# Make and save visuals
p.clusters <- SCpubr::do_DimPlot(rau.sn,
                   label = T,
                   repel = T,
                   font.size = 14,
                   label.size = 4) + 
                 NoLegend() + 
                 ggtitle("Clusters at 0.2 resolution with RunHarmony") +
                 theme(
                   plot.title = element_text(size = 14, hjust = 0.5))
p.origin <- SCpubr::do_DimPlot(rau.sn,
                               group.by = "orig.ident",
                               label = T,
                               repel = T,
                               font.size = 14,
                               label.size = 4) + 
                      NoLegend() + 
                      ggtitle("Distribution of samples") +
                      theme(
                        plot.title = element_text(size = 14, hjust = 0.5))

# Function to plot features of clusters as ridge density plots
makeRidgePlot <- function(feature, title){
rau.sn@meta.data |> 
  ggplot(aes(x = get(feature), y = seurat_clusters, fill = after_stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = as.character(feature),
                     option = "viridis",
                     direction = -1
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(y = "Seurat Cluster",
       x = title)
  }

# Apply the function
p.features <- makeRidgePlot("nFeature_RNA", "Feature count of droplets")
p.count <- makeRidgePlot("nCount_RNA", "Read count within droplets")
p.mito <- makeRidgePlot("PercentMito", "Percent of mitocondrial transcripts")
p.doublet <- makeRidgePlot("DoubletScore", "Doublet likelihood score")


# Wrap the plots

design <- "EAB
           FCD"


png(file = "data/processed/single_cell/merges/rau_patterson/09132023/cluster_features_3.png",
    width = 1920, 
    height = 1080,
    units = "px")

plots <- patchwork::wrap_plots(A = p.features, B = p.count,
                               C = p.mito,     D = p.doublet, 
                               E = p.clusters, F = p.origin, 
                               design = design)

plots

dev.off()

