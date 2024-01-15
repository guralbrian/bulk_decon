---
title: "Rau Doublet/Ambient RNA Removal"
author: "Brian Gural"
date: "2023-08-14"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


```{r load libs, message=FALSE, echo = T, warning=FALSE, results= 'hide', cache=FALSE, include=F}
# List libraries
setwd("/proj/raulab/users/brian/r_projects/cc_composition")
libs <- c("Seurat", "ggplot2", "DESeq2", "patchwork","SeuratDisk", "reshape2",
          "tidyverse", "SingleCellExperiment", "harmony", "SCpubr", 
          "AUCell", "viridis", "gplots", "scales", "ggrepel", "gridExtra", "scCustomize",
          "httr","matrixStats", "scran", "scuttle", "scater", "DropletUtils", "scDblFinder") # list libraries here
# Require all of them
lapply(libs, require, character.only = T)
rm(libs)
# Also load in-house functions
source("jensen/scripts/functions/decon_all.R")
```


```{r load data, message=FALSE, echo = T, warning=FALSE, cache=FALSE, include=F}
# Load Seurats as separate batches
sn.1 <- LoadH5Seurat("jensen/data/raw/single_cell/patterson/seurat/b6_1.h5seurat")

sn.1 <- subset(sn.1, nCount_RNA > 0)
sn.1$UMI <- colnames(sn.1)

sn.2 <- LoadH5Seurat("jensen/data/raw/single_cell/patterson/seurat/b6_2.h5seurat")
```
## Rational

This script is meant to remove single cell doublets and ambient RNA from a snRNAseq dataset produced by Christoph Rau. Each sample is the aggregate of 5 samples, 

```{r subset empty drops}

empty.drops <- sn.1 |>
               subset(nCount_RNA <= 300)

empty.subset <- sample(colnames(empty.drops), 10000)


keep.drops <- sn.1 |>
              subset(nCount_RNA >= 300 | UMI %in% empty.subset)

```

```{r cars}
  # Find empty droplets 
  # make sce
  sce.1 <- as.SingleCellExperiment(keep.drops)
  
# rank each cell by the # counts
  bcrank <- barcodeRanks(counts(sce.1))
  uniq <- !duplicated(bcrank$rank)
  
# Now plot the ranks with the counts, look for the knee and inflection points
  plot(bcrank$rank[uniq], bcrank$total[uniq], log="y", xlab="Rank", ylab="Total UMI count", cex.lab=1.2) 
  o <- order(bcrank$rank)
    lines(bcrank$rank[o], bcrank$fitted[o], col="red")
    abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2) 
    abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2) 
    abline(h=500, col="darkorange", lty=2) 
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

```


```{r emptyDrops limited}

  limit <- metadata(bcrank)$inflection   
  e.out <- emptyDrops(counts(sce.1), test.ambient=TRUE, lower=500)
  
  # Large p-value peak near 0, some cells are likely being considered as ambient RNA
  # Rerun with limit cuttoff to exclude them
  hist(e.out$PValue[e.out$Total <= limit & e.out$Total > 0],
    xlab="P-value", main="", col="grey80") 

```

```{r QC}
# limit to droplets lower than the false discovery threshold 

sce.1.2 <- sce.1[,which(e.out$FDR <= 0.01)]
  
library(scuttle)

# find mitocondrial gene percent outliers and exclude
is.mito <- grep("^mt-", rownames(sce.1.2))
qc <- perCellQCMetrics(sce.1.2, subsets=list(MT=is.mito))
discard.mito <- isOutlier(qc$subsets_MT_percent, type="higher")
summary(discard.mito)


# Plot 
plot(qc$sum, qc$subsets_MT_percent, log="xy",
    xlab="Total count", ylab='Mitochondrial %')
abline(h=attr(discard.mito, "thresholds")["higher"], col="red")


colData(sce.1.2) <- cbind(colData(sce.1.2), qc)
sce.1.2 <- sce.1.2[,!discard.mito]
```


```{r QC cluster/variance}
  
  # Single cell processing (clustering and normalization)
  # References http://bioconductor.org/books/3.14/OSCA.workflows/unfiltered-human-pbmcs-10x-genomics.html#unfiltered-human-pbmcs-10x-genomics
  clusters <- quickCluster(sce.1.2)
  sce.1.2 <- computeSumFactors(sce.1.2, cluster=clusters)
  sce.1.2 <- logNormCounts(sce.1.2)
  
  
  # Variance modeling
  sce.dec <- modelGeneVarByPoisson(sce.1.2)
  sec.top <- getTopHVGs(sce.dec, prop=0.1)
  
  # Per-gene variance as a function of the mean for the log-expression values in the PBMC dataset. 
  # Each point represents a gene (black) with the mean-variance trend (blue) fitted to simulated Poisson counts.
  plot(sce.dec$mean, sce.dec$total, pch=16, cex=0.5,
    xlab="Mean of log-expression", ylab="Variance of log-expression")
  curfit <- metadata(sce.dec)
  curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
```

```{r QC dimensional reduction }
  
set.seed(10000)
sce.1.2 <- denoisePCA(sce.1.2, subset.row=sec.top, technical=sce.dec)

set.seed(100000)
sce.1.2 <- runTSNE(sce.1.2, dimred="PCA")

set.seed(1000000)
sce.1.2 <- runUMAP(sce.1.2, dimred="PCA")

# Not much variance is being explained by these PCs...

plot(attr(reducedDim(sce.1.2), "varExplained"), xlab="PC", ylab="Variance explained (%)")

```


```{r QC clustering}

g <- buildSNNGraph(sce.1.2, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.1.2) <- factor(clust)


plotTSNE(sce.1.2, colour_by="label")
```

```{r remove ambience}

amb <- metadata(e.out)$ambient[,1]
stripped <- sce.1.2[names(amb),]

out <- removeAmbience(counts(stripped), ambient=amb, groups=colLabels(stripped))
```

```{r remove ambience visual}
library(scater)
counts(stripped, withDimnames=FALSE) <- out
stripped <- logNormCounts(stripped)
plotAmbience <- function(gene){
gridExtra::grid.arrange(
            plotExpression(sce.1.2, x="label", colour_by="label", features=gene) + 
        ggtitle("Before"),
          plotExpression(stripped, x="label", colour_by="label", features=gene) + 
        ggtitle("After"),ncol=2) }

genes <- c("Col1a1","Acta2", "Myh6", "Tnnt2", "D830005E20Rik") 

lapply(genes, function(x){
  plotAmbience(x)
})
```

```{r doublets and visuals}
  # Find Doublets ####
  
  dbl.dens <- computeDoubletDensity(stripped, #subset.row=top.mam, 
                                    d=ncol(reducedDim(stripped)),subset.row=sec.top)
  
  stripped$DoubletScore <- dbl.dens
  
  # Save plotUMAP with DoubletScore
  
  scater::plotUMAP(stripped, colour_by="DoubletScore")
  # Save plotColData with DoubletScore
  plotColData(stripped, x="label", y="DoubletScore", colour_by="label")+
    geom_hline(yintercept = quantile(colData(stripped)$DoubletScore, 0.95),lty="dashed",color="red")
  
  cut_off <- quantile(stripped$DoubletScore,0.95)
  stripped$isDoublet <- c("no","yes")[factor(as.integer(stripped$DoubletScore>=cut_off),levels=c(0,1))]
  
  sn.clean <- as.Seurat(stripped)
  SaveH5Seurat(sn.clean, "jensen/results/doublets_ambient_rna/08142023/B6_1")
```


```{r load data, message=FALSE, echo = T, warning=FALSE, cache=FALSE, include=F}
sn.2 <- LoadH5Seurat("jensen/data/raw/single_cell/patterson/seurat/b6_2.h5seurat")

sn.2 <- subset(sn.2, nCount_RNA > 0)
sn.2$UMI <- colnames(sn.2)
```
## Rational

This script is meant to remove single cell doublets and ambient RNA from a snRNAseq dataset produced by Christoph Rau. Each sample is the aggregate of 5 samples, 

```{r subset empty drops}

empty.drops <- sn.2 |>
               subset(nCount_RNA <= 300)

empty.subset <- sample(colnames(empty.drops), 10000)


keep.drops <- sn.2 |>
              subset(nCount_RNA >= 300 | UMI %in% empty.subset)

```

```{r cars}
  # Find empty droplets 
  # make sce
  sce.2 <- as.SingleCellExperiment(keep.drops)
  
# rank each cell by the # counts
  bcrank <- barcodeRanks(counts(sce.2))
  uniq <- !duplicated(bcrank$rank)
  
# Now plot the ranks with the counts, look for the knee and inflection points
  plot(bcrank$rank[uniq], bcrank$total[uniq], log="y", xlab="Rank", ylab="Total UMI count", cex.lab=1.2) 
  o <- order(bcrank$rank)
    lines(bcrank$rank[o], bcrank$fitted[o], col="red")
    abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2) 
    abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2) 
    abline(h=500, col="darkorange", lty=2) 
    legend("bottomleft", legend=c("Inflection", "Knee"), col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
  
```


```{r emptyDrops limited}

  limit <- metadata(bcrank)$inflection   
  e.out <- emptyDrops(counts(sce.2), test.ambient=TRUE, lower=500)
  
  # Large p-value peak near 0, some cells are likely being considered as ambient RNA
  # Rerun with limit cuttoff to exclude them
  hist(e.out$PValue[e.out$Total <= limit & e.out$Total > 0],
    xlab="P-value", main="", col="grey80") 

```

```{r QC}
# limit to droplets lower than the false discovery threshold 

sce.2.2 <- sce.2[,which(e.out$FDR <= 0.01)]
  
library(scuttle)

# find mitocondrial gene percent outliers and exclude
is.mito <- grep("^mt-", rownames(sce.2.2))
qc <- perCellQCMetrics(sce.2.2, subsets=list(MT=is.mito))
discard.mito <- isOutlier(qc$subsets_MT_percent, type="higher")
summary(discard.mito)


# Plot 
plot(qc$sum, qc$subsets_MT_percent, log="xy",
    xlab="Total count", ylab='Mitochondrial %')
abline(h=attr(discard.mito, "thresholds")["higher"], col="red")


colData(sce.2.2) <- cbind(colData(sce.2.2), qc)
sce.2.2 <- sce.2.2[,!discard.mito]
```


```{r QC cluster/variance}
  
  # Single cell processing (clustering and normalization)
  # References http://bioconductor.org/books/3.14/OSCA.workflows/unfiltered-human-pbmcs-10x-genomics.html#unfiltered-human-pbmcs-10x-genomics
  clusters <- quickCluster(sce.2.2)
  sce.2.2 <- computeSumFactors(sce.2.2, cluster=clusters)
  sce.2.2 <- logNormCounts(sce.2.2)
  
  
  # Variance modeling
  sce.dec <- modelGeneVarByPoisson(sce.2.2)
  sec.top <- getTopHVGs(sce.dec, prop=0.1)
  
  # Per-gene variance as a function of the mean for the log-expression values in the PBMC dataset. 
  # Each point represents a gene (black) with the mean-variance trend (blue) fitted to simulated Poisson counts.
  plot(sce.dec$mean, sce.dec$total, pch=16, cex=0.5,
    xlab="Mean of log-expression", ylab="Variance of log-expression")
  curfit <- metadata(sce.dec)
  curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
```

```{r QC dimensional reduction }
  
set.seed(10000)
sce.2.2 <- denoisePCA(sce.2.2, subset.row=sec.top, technical=sce.dec)

set.seed(100000)
sce.2.2 <- runTSNE(sce.2.2, dimred="PCA")

set.seed(1000000)
sce.2.2 <- runUMAP(sce.2.2, dimred="PCA")
```


```{r QC clustering}

g <- buildSNNGraph(sce.2.2, k=20, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.2.2) <- factor(clust)


plotTSNE(sce.2.2, colour_by="label")
```

```{r remove ambience}

amb <- metadata(e.out)$ambient[,1]
stripped <- sce.2.2[names(amb),]

out.2 <- removeAmbience(counts(stripped), ambient=amb, groups=colLabels(stripped))
```

```{r remove ambience visual}
library(scater)
counts(stripped, withDimnames=FALSE) <- out.2
stripped <- logNormCounts(stripped)
plotAmbience <- function(gene){
gridExtra::grid.arrange(
            plotExpression(sce.2.2, x="label", colour_by="label", features=gene) + 
        ggtitle("Before"),
          plotExpression(stripped, x="label", colour_by="label", features=gene) + 
        ggtitle("After"),ncol=2) }

genes <- c("Col1a1","Acta2", "Myh6", "Tnnt2", "D830005E20Rik") 

lapply(genes, function(x){
  plotAmbience(x)
})
```













```{r doublets and visuals}
  # Find Doublets ####
  
  dbl.dens <- computeDoubletDensity(stripped, #subset.row=top.mam, 
                                    d=ncol(reducedDim(stripped)),subset.row=sec.top)
  
  stripped$DoubletScore <- dbl.dens
  
  # Save plotUMAP with DoubletScore
  
  scater::plotUMAP(stripped, colour_by="DoubletScore")
  # Save plotColData with DoubletScore
  plotColData(stripped, x="label", y="DoubletScore", colour_by="label")+
    geom_hline(yintercept = quantile(colData(stripped)$DoubletScore, 0.95),lty="dashed",color="red")
  
  cut_off <- quantile(stripped$DoubletScore,0.95)
  stripped$isDoublet <- c("no","yes")[factor(as.integer(stripped$DoubletScore>=cut_off),levels=c(0,1))]
  
  sn.clean <- as.Seurat(stripped)
  SaveH5Seurat(sn.clean, "jensen/results/doublets_ambient_rna/08142023/B6_2")
```