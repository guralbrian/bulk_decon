
# List libraries
setwd("/proj/raulab/users/brian/r_projects/bulk_decon")
libs <- c("Seurat", "ggplot2", "DESeq2", "patchwork","SeuratDisk", "reshape2",
          "tidyverse", "SingleCellExperiment", "harmony", "SCpubr", 
          "AUCell", "viridis", "gplots", "scales", "ggrepel", "gridExtra", "scCustomize",
          "httr","matrixStats", "scran", "scuttle", "scater", "DropletUtils", "scDblFinder") # list libraries here
# Require all of them
lapply(libs, require, character.only = T)

# Also load in-house functions
source("jensen/scripts/functions/decon_all.R")


sn.1<- LoadH5Seurat("jensen/data/raw/single_cell/patterson/seurat/b6_1.h5seurat")

sn.2 <- LoadH5Seurat("jensen/data/raw/single_cell/patterson/seurat/b6_2.h5seurat")


# Find empty droplets 
# make sce
sce.1 <- as.SingleCellExperiment(sn.1)

bcrank <- barcodeRanks(counts(sce.1))
uniq <- !duplicated(bcrank$rank)

# Now at each plot, you wrap the plot call with png(), and dev.off() to save the image
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Save plot
ggsave(filename = "jensen/results/doublets_ambient_rna/08142023/B6_1_elbow.png", units = "px",dpi=300, width = 1000, height = 900 )

# set cutoff limit and consider droplets below it to be empty
set.seed(100)
e.out <- emptyDrops(counts(sce.1), test.ambient=TRUE)

# Large p-value peak near 0, some cells are likely being considered as ambient RNA
# Rerun with limit cuttoff to exclude them
hist(e.out$PValue[ e.out$Total > 0],
     xlab="P-value", main="", col="grey80") 



limit <- metadata(bcrank)$inflection   
e.out <- emptyDrops(counts(sce.1), test.ambient=TRUE, lower=limit)

# Large p-value peak near 0, some cells are likely being considered as ambient RNA
# Rerun with limit cuttoff to exclude them

sce.1.2 <- sce.1[,which(e.out$FDR <= 0.01)]


# Define and remove ambient RNA 
clusters <- quickCluster(sce.1.2)
sce.1.2 <- computeSumFactors(sce.1.2, cluster=clusters)

# evaluate ambinat RNA contamination in the empty droplets
amb <- metadata(e.out)$ambient[,1]

sce.1.2 <- logNormCounts(sce.1.2)
set.seed(1000)
# modeling variables
dec.pbmc <- modelGeneVarByPoisson(sce.1.2)
# calcualte top features
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)
#
set.seed(1000)
# Evaluate PCs
sce.1.2 <- denoisePCA(sce.1.2, subset.row=top.pbmc, technical=dec.pbmc)
# make UMAP plot
sce.1.2 <- runUMAP(sce.1.2, dimred="PCA")
g <- buildSNNGraph(sce.1.2, k=25, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.1.2) <- factor(clust)

# Save plotUMAP

scater::plotUMAP(sce.1.2, colour_by="label")
ggsave(filename = "jensen/results/doublets_ambient_rna/08142023/B6_1_umap.png", units = "px", dpi =300, width = 1000, height = 900)


stripped <- sce.1.2[names(amb),]
out <- removeAmbience(counts(stripped), ambient= amb, groups = colLabels(stripped))
# load correccted counts into scce object
counts(stripped, withDimnames=FALSE) <- out
stripped <- logNormCounts(stripped)


# Find Doublets ####

dbl.dens <- computeDoubletDensity(stripped, #subset.row=top.mam, 
                                  d=ncol(reducedDim(stripped)),subset.row=top.pbmc)

stripped$DoubletScore <- dbl.dens

# Save plotUMAP with DoubletScore

scater::plotUMAP(stripped, colour_by="DoubletScore")
ggsave(filename = "jensen/results/doublets_ambient_rna/08142023/B6_1_doub_umap.png", units = "px",dpi=300, width = 1000, height = 900 )

# Save plotColData with DoubletScore
plotColData(stripped, x="label", y="DoubletScore", colour_by="label")+
  geom_hline(yintercept = quantile(colData(stripped)$DoubletScore, 0.95),lty="dashed",color="red")
ggsave(filename = "jensen/results/doublets_ambient_rna/08142023/B6_1_doub_score.png", units = "px",dpi=300, width = 1000, height = 900 )


cut_off <- quantile(stripped$DoubletScore,0.95)
stripped$isDoublet <- c("no","yes")[factor(as.integer(stripped$DoubletScore>=cut_off),levels=c(0,1))]

sn.clean <- as.Seurat(stripped)
SaveH5Seurat(sn.clean, "jensen/results/doublets_ambient_rna/08142023/B6_1")


# Find empty droplets 
# make sce
sce.1 <- as.SingleCellExperiment(sn.2)

bcrank <- barcodeRanks(counts(sce.1))
uniq <- !duplicated(bcrank$rank)

# Now at each plot, you wrap the plot call with png(), and dev.off() to save the image
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Save plot
ggsave(filename = "jensen/results/doublets_ambient_rna/08142023/B6_1_elbow.png", units = "px",dpi=300, width = 1000, height = 900 )

# set cutoff limit and consider droplets below it to be empty
set.seed(100)
e.out <- emptyDrops(counts(sce.1), test.ambient=TRUE)

# Large p-value peak near 0, some cells are likely being considered as ambient RNA
# Rerun with limit cuttoff to exclude them
hist(e.out$PValue[ e.out$Total > 0],
     xlab="P-value", main="", col="grey80") 



limit <- metadata(bcrank)$inflection   
e.out <- emptyDrops(counts(sce.1), test.ambient=TRUE, lower=limit)

# Large p-value peak near 0, some cells are likely being considered as ambient RNA
# Rerun with limit cuttoff to exclude them
hist(e.out$PValue[e.out$Total <= limit & e.out$Total > 0],
     xlab="P-value", main="", col="grey80") 



sce.1.2 <- sce.1[,which(e.out$FDR <= 0.01)]


# Define and remove ambient RNA 
clusters <- quickCluster(sce.1.2)
sce.1.2 <- computeSumFactors(sce.1.2, cluster=clusters)

# evaluate ambinat RNA contamination in the empty droplets
amb <- metadata(e.out)$ambient[,1]

sce.1.2 <- logNormCounts(sce.1.2)
set.seed(1000)
# modeling variables
dec.pbmc <- modelGeneVarByPoisson(sce.1.2)
# calcualte top features
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)
#
set.seed(1000)
# Evaluate PCs
sce.1.2 <- denoisePCA(sce.1.2, subset.row=top.pbmc, technical=dec.pbmc)
# make UMAP plot
sce.1.2 <- runUMAP(sce.1.2, dimred="PCA")
g <- buildSNNGraph(sce.1.2, k=25, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.1.2) <- factor(clust)

# Save plotUMAP

scater::plotUMAP(sce.1.2, colour_by="label")
ggsave(filename = "jensen/results/doublets_ambient_rna/08142023/B6_1_umap.png", units = "px", dpi =300, width = 1000, height = 900)


stripped <- sce.1.2[names(amb),]
out <- removeAmbience(counts(stripped), ambient= amb, groups = colLabels(stripped))
# load correccted counts into scce object
counts(stripped, withDimnames=FALSE) <- out
stripped <- logNormCounts(stripped)

# Find Doublets ####

dbl.dens <- computeDoubletDensity(stripped, #subset.row=top.mam, 
                                  d=ncol(reducedDim(stripped)),subset.row=top.pbmc)

stripped$DoubletScore <- dbl.dens

# Save plotUMAP with DoubletScore

scater::plotUMAP(stripped, colour_by="DoubletScore")
ggsave(filename = "jensen/results/doublets_ambient_rna/08142023/B6_1_doub_umap.png", units = "px",dpi=300, width = 1000, height = 900 )

# Save plotColData with DoubletScore
plotColData(stripped, x="label", y="DoubletScore", colour_by="label")+
  geom_hline(yintercept = quantile(colData(stripped)$DoubletScore, 0.95),lty="dashed",color="red")
ggsave(filename = "jensen/results/doublets_ambient_rna/08142023/B6_1_doub_score.png", units = "px",dpi=300, width = 1000, height = 900 )


cut_off <- quantile(stripped$DoubletScore,0.95)
stripped$isDoublet <- c("no","yes")[factor(as.integer(stripped$DoubletScore>=cut_off),levels=c(0,1))]

sn.clean <- as.Seurat(stripped)
SaveH5Seurat(sn.clean, "jensen/results/doublets_ambient_rna/08142023/B6_2")
