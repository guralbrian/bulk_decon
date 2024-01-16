


libs <- c("Seurat", "ggplot2", "DESeq2", "patchwork","SeuratDisk", "reshape2",
          "tidyverse", "SingleCellExperiment", "harmony", "SCpubr", 
          "AUCell", "viridis", "gplots", "scales", "ggrepel", "gridExtra", "scCustomize",
          "httr","matrixStats", "scran", "scuttle", "scater", "DropletUtils", "scDblFinder") # list libraries here

libs <- c("Seurat", "SeuratDisk", "scran", "DropletUtils", "scater", "scDblFinder")

# Require all of them
lapply(libs, require, character.only = T)


sample_name <- "b6_2"
# Load Seurat as separate batches
sn.1 <- LoadH5Seurat(paste0("data/processed/single_cell/", sample_name, ".h5seurat"))



#### Filter droplets ####
sn.1 <- subset(sn.1, nCount_RNA > 0)
sn.1$UMI <- colnames(sn.1)
# Find low UMI droplets (putative empty drops) and include a subset of them
empty.drops <- names(sn.1$nCount_RNA[which(sn.1$nCount_RNA <= 300)])
sn.1 <- sn.1 |>
        subset(nCount_RNA >= 300 | UMI %in% sample(empty.drops, 10000))

# Convert to sce
sce.1 <- as.SingleCellExperiment(sn.1)
  
# rank each cell by the # counts
bcrank <- barcodeRanks(counts(sce.1))
uniq <- !duplicated(bcrank$rank)
  
# Save plot
png_path <- paste0("results/2_ambient_doublets/",sample_name, "_elbow.png")
png(png_path, width = 14, height = 8, units = "in", res = 300)

# Plot the ranks with the counts, look for the knee and inflection points
plot(bcrank$rank[uniq], bcrank$total[uniq], log="y", xlab="Rank", ylab="Total UMI count", cex.lab=1.2) 
o <- order(bcrank$rank)
lines(bcrank$rank[o], bcrank$fitted[o], col="red")
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2) 
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2) 
abline(h=500, col="darkorange", lty=2) 
legend("bottomleft", legend=c("Inflection", "Knee"), col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
    
dev.off()  
    
# set cutoff limit and consider droplets below it to be empty
# Large p-value peak near 0, some cells are likely being considered as ambient RNA
# Rerun with limit cuttoff to exclude them
set.seed(100)
limit <- metadata(bcrank)$inflection   
e.out <- emptyDrops(counts(sce.1), test.ambient=TRUE, lower=500)

# limit to droplets lower than the false discovery threshold 
sce.1 <- sce.1[,which(e.out$FDR <= 0.01)]

# find mitocondrial gene percent outliers and exclude
is.mito <- grep("^mt-", rownames(sce.1))
qc <- perCellQCMetrics(sce.1, subsets=list(MT=is.mito))
discard.mito <- isOutlier(qc$subsets_MT_percent, type="higher")

# Save plot
png_path <- paste0("results/2_ambient_doublets/",sample_name, "_mito.png")
png(png_path, width = 14, height = 8, units = "in", res = 300)

# Plot the mitocontrial distribution
plot(qc$sum, qc$subsets_MT_percent, log="xy",
     xlab="Total count", ylab='Mitochondrial %')
abline(h=attr(discard.mito, "thresholds")["higher"], col="red")

dev.off()  

# Add in this meta data, remove the putative mito contaminated droplets
colData(sce.1) <- cbind(colData(sce.1), qc)
sce.1 <- sce.1[,!discard.mito]

#### Single cell processing (clustering and normalization) ####
# References http://bioconductor.org/books/3.14/OSCA.workflows/unfiltered-human-pbmcs-10x-genomics.html#unfiltered-human-pbmcs-10x-genomics
clusters <- quickCluster(sce.1)
sce.1 <- computeSumFactors(sce.1, cluster=clusters)
sce.1 <- logNormCounts(sce.1)
  
# Variance modeling
sce.dec <- modelGeneVarByPoisson(sce.1)
sec.top <- getTopHVGs(sce.dec, prop=0.1)
  
# Per-gene variance as a function of the mean for the log-expression values
# Each point represents a gene (black) with the mean-variance trend (blue) fitted to simulated Poisson counts.
  # Save plot
  png_path <- paste0("results/2_ambient_doublets/",sample_name, "_gene_var.png")
  png(png_path, width = 14, height = 8, units = "in", res = 300)
  
  plot(sce.dec$mean, sce.dec$total, pch=16, cex=0.5,
       xlab="Mean of log-expression", ylab="Variance of log-expression")
  curfit <- metadata(sce.dec)
  curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
  
  dev.off()  

# Reduce dimensionality and project
set.seed(100)
sce.1 <- denoisePCA(sce.1, subset.row = sec.top, technical = sce.dec) |>
  scater::runTSNE(dimred="PCA") |> 
  scater::runUMAP(dimred="PCA")

# Plot variance explained by each PC
png_path <- paste0("results/2_ambient_doublets/",sample_name, "_pca_var.png")
png(png_path, width = 14, height = 8, units = "in", res = 300)

plot(attr(reducedDim(sce.1), "varExplained"), xlab="PC", ylab="Variance explained (%)")

dev.off()  

# Cluster the reduced droplet data
g <- buildSNNGraph(sce.1, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.1) <- factor(clust)

png_path <- paste0("results/2_ambient_doublets/",sample_name, "_tnse.png")
png(png_path, width = 14, height = 8, units = "in", res = 300)

# plot the clusters
plotTSNE(sce.1, colour_by="label")

dev.off()  


#### Remove ambient RNA ####
amb <- metadata(e.out)$ambient[,1]
stripped <- sce.1[names(amb),]
out <- removeAmbience(counts(stripped), ambient=amb, groups=colLabels(stripped))

# Visualize the before/after changes
counts(stripped, withDimnames=FALSE) <- out
stripped <- logNormCounts(stripped)
plotAmbience <- function(gene){
  

gridExtra::grid.arrange(
            plotExpression(sce.1, x="label", colour_by="label", features=gene) + 
        ggtitle("Before"),
          plotExpression(stripped, x="label", colour_by="label", features=gene) + 
        ggtitle("After"),ncol=2) }

# These are example genes. They should be highly expressed only in cardiomyocytes
genes <- c("Col1a1","Acta2", "Myh6", "Tnnt2") 

lapply(genes, function(x){
  png_path <- paste0("results/2_ambient_doublets/", sample_name, "_", x, "_ambient.png")
  png(png_path, width = 14, height = 8, units = "in", res = 300)
  plotAmbience(x)
  dev.off()  
})


#### Find Doublets ####
dbl.dens <- scDblFinder::computeDoubletDensity(stripped, #subset.row=top.mam, 
                d=ncol(reducedDim(stripped)),subset.row=sec.top)
  
stripped$DoubletScore <- dbl.dens

# Save plotUMAP with DoubletScore

png_path <- paste0("results/2_ambient_doublets/", sample_name, "_doublets.png")
png(png_path, width = 14, height = 8, units = "in", res = 300)
gridExtra::grid.arrange(
  scater::plotUMAP(stripped, colour_by="DoubletScore"),
  plotColData(stripped, x="label", y="DoubletScore", colour_by="label")+
    geom_hline(yintercept = quantile(colData(stripped)$DoubletScore, 0.95),lty="dashed",color="red"),
  ncol = 2)
dev.off()  

# Call doublets by quantile
cut_off <- quantile(stripped$DoubletScore,0.95)
stripped$isDoublet <- c("no","yes")[factor(as.integer(stripped$DoubletScore>=cut_off),levels=c(0,1))]

# Save data
SaveH5Seurat(as.Seurat(stripped), 
             paste0("data/processed/single_cell/", sample_name, "_no_doublets"))
