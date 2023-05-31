# tabula muris
install.packages("dplyr")
BiocManager::install("TabulaMurisSenisData")
library("TabulaMurisSenisData")
listTabulaMurisSenisTissues("Droplet")

sn_muris <- TabulaMurisSenisDroplet(
  tissues = "Heart_and_Aorta",
  processedCounts = FALSE,
  reducedDims = TRUE,
  infoOnly = FALSE
)

sn_muris <- sn_muris[[1]]

# Realize the DelayedMatrix
sn_muris_matrix <- as.matrix(assay(sn_muris, "counts"))

# Create a Seurat object
sn_muris_seurat <- Seurat::CreateSeuratObject(counts = sn_muris_matrix)
# Add metadata from SingleCellExperiment to Seurat

for(i in colnames(colData(sn_muris))){
sn_muris_seurat <- Seurat::AddMetaData(object = sn_muris_seurat, col.name = i, metadata = colData(sn_muris)[[i]])
}


feat.ft <- c("PercentMito", "nCount_RNA", "PercentRibo")
vln.ft <- c("PercentMito")

sn_muris_seurat[["PercentMito"]] <- PercentageFeatureSet(object = sn_muris_seurat, pattern = "^Mt-")
sn_muris_seurat$PercentRibo      <- 0
sn_muris_seurat$origin <- "tabula_muris"
sn.clust$origin <-  "rau"

sn.all <- merge(sn, sn_muris_seurat)
sn.all$origin <- as.factor(sn.all$origin)


subset = T
min.rna.ft = 200
max.rna.ft = 2500
min.rna.ct = 800
max.mt.pt  = 0.05
regress.by = "orig.ident"
harmony    = T
res        = 0.2
nfeatures  = 2000

max.rna.ft = 4000
max.mt.pt = 5 
res = 0.1 
regress.by = "origin"


sn.clust.1 <- subset(sn.all, 
                     subset = nFeature_RNA   > min.rna.ft     & 
                       nFeature_RNA   < max.rna.ft     &
                       nCount_RNA     > min.rna.ct     &
                       PercentMito   <= max.mt.pt)

seurat.obj <- sn.clust.1 |>
  NormalizeData(verbose = F) |>
  FindVariableFeatures(verbose = F, nfeatures = nfeatures) |>
  ScaleData(verbose = F) |>
  RunPCA(verbose = F) 

if(harmony == T){
  seurat.obj@meta.data[[regress.by]] <- droplevels(seurat.obj@meta.data[[regress.by]]) 
  seurat.obj <- RunHarmony(seurat.obj, 
                           group.by.vars = regress.by,
                           verbose = F)
}

# find elbow
# Determine percent of variation associated with each PC
pct <- seurat.obj[["pca"]]@stdev / sum(seurat.obj[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

if(harmony == T){
  seurat.obj <- seurat.obj |>
    FindNeighbors(dims = 1:pcs, reduction = "harmony", verbose = F) |>
    FindClusters(resolution = res, verbose = F) |>
    RunUMAP(dims = 1:pcs, reduction = "harmony", verbose = F)
}


seurat.obj <- seurat.obj |>
  FindNeighbors(dims = 1:pcs, reduction = "harmony", verbose = F) |>
  FindClusters(resolution = res, verbose = F) |>
  RunUMAP(dims = 1:pcs, reduction = "harmony", verbose = F)

sn.clust.1 <- seurat.obj

sn.clust.2 <- sn.clust.1

# music gene weights
gene_weights <- music_basis(sn.sce,clusters = "ident", samples = "orig.ident")

# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(gene_weights$Disgn.mtx + 1e-6)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d <- dist(t(log(gene_weights$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
# hc2 <- hclust(d, method = "complete" )
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')
gene_sigma <- gene_weights$Sigma


# make MitoClusters() only consider nuclei from Rau

data <- sn.clust.1
cutoff = 3
origins.considered = "rau"
  sn.mito <- subset(data, subset = PercentMito > cutoff &
                                   origin %in% origins.considered)
  sn.mito.tb <- table(sn.mito$seurat_clusters)
  sn.origin <- subset(data, subset = origin %in% origins.considered)
  sn.tb <- table(sn.origin$seurat_clusters)
  clust.pcs <- (sn.mito.tb / sn.tb) * 100 |>
    round(digits = 1)
  return(clust.pcs)
  
  
  
  # troubleshoot ClusterSeurat ####
  
  seurat_obj <- sn_all
  subset = T
  min.rna.ft = 200
  max.rna.ft = 2500
  min.rna.ct = 800
  max.mt.pt  = 0.05
  regress.by = "orig.ident"
  harmony    = T
  res        = 0.2
  nfeatures  = 2000
  max.rna.ft = 4000
  max.mt.pt = 1
  res = 0.05
  regress.by = c("origin", "PercentMito")
  
  
  seurat_obj <- subset(seurat_obj, 
                           subset = nFeature_RNA   > min.rna.ft     & 
                             nFeature_RNA   < max.rna.ft     &
                             nCount_RNA     > min.rna.ct     &
                             PercentMito   <= max.mt.pt)

    seurat_obj <- seurat_obj |>
      NormalizeData(verbose = F) |>
      FindVariableFeatures(verbose = F, nfeatures = nfeatures) |>
      ScaleData(verbose = F) |>
      RunPCA(verbose = F) 
    
    if(harmony == T){
      seurat_obj@meta.data[[regress.by]] <- droplevels(seurat_obj@meta.data[[regress.by]]) 
      seurat_obj <- RunHarmony(seurat_obj, 
                               group.by.vars = "orig.ident",
                               verbose = F,
                               project.dim = F)
    }
    
    # find elbow
    # Determine percent of variation associated with each PC
    pct <- seurat.obj[["pca"]]@stdev / sum(seurat.obj[["pca"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    pcs <- min(co1, co2)

      seurat.obj <- seurat.obj |>
        FindNeighbors(dims = 1:pcs, reduction = "harmony", verbose = F) |>
        FindClusters(resolution = res, verbose = F) |>
        RunUMAP(dims = 1:pcs, reduction = "harmony", verbose = F)
      