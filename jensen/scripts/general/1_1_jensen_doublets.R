# 1_1_jensen_qc

# this script is meant to improve the pre-processing of the individual single cell datasets used in the pipeline
# This specifc part is meant to load each of the seurats loading in 1_jensen_load.R then filter/scale/and normalize each
# This will also implement two new methods: DoubletFinder and SoupX to remove doublets and ambient RNA contamination
# each of these will be run seperatel and saved in their own directory together 


# Load data


# load libraries
libs <- c("Seurat", "ggplot2", "DESeq2", "patchwork","SeuratDisk", "MuSiC", "reshape2",
          "tidyverse", "SingleCellExperiment","harmony", "SCpubr","shiny", 
          "AUCell", "viridis", "gplots", "scales", "ggrepel", "gridExtra", "scCustomize",
          "httr","readxl","matrixStats", "TabulaMurisSenisData") # list libraries here
lapply(libs, require, character.only = T)
source("jensen/scripts/functions/decon_all.R")

# Single cell
sn.rau <- LoadH5Seurat("jensen/data/raw/b6_sn.h5seurat")
sn.tm  <- LoadH5Seurat("jensen/data/processed/tabula_muris.h5seurat")
sn.tm$orig.ident <- paste0("tm_", sn.tm$mouse.id)
sn.tm <- subset(sn.tm, tissue_free_annotation == "Heart")

sn.wu  <- LoadH5Seurat("jensen/data/processed/single_cell/wu_2021.h5seurat")
sn.mt  <- LoadH5Seurat("jensen/data/processed/single_cell/martini_2019.h5seurat")

sn.rau$origin <- "rau"
sn.mt$origin <- "martini"


sn.rau$PercentMito <- PercentageFeatureSet(sn.rau, pattern = "^mt-")
sn.wu$PercentMito <- PercentageFeatureSet(sn.wu, pattern = "^mt-")
sn.mt$PercentMito <- PercentageFeatureSet(sn.mt, pattern = "^mt-")

# QC on each

# New ClusterSeurat includes SoupX and DoubletFinder
# and considers distribution of counts and features of each origin
sub.rau <- unique(sn.rau$orig.ident)

sn.rau <- lappy(uniqueClusterByMeta(sn.rau, sub.rau))

# Save in new fotlder
table(sn.rau$DF.classifications_0.25_0.07_1747)



# finding empty droplets
library("DropletUtils")

# make sce

sce.rau <- as.SingleCellExperiment(sn.rau)

bcrank <- barcodeRanks(counts(sce.rau))
uniq <- !duplicated(bcrank$rank)

plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
set.seed(100)
limit <- 300   
e.out <- emptyDrops(counts(sce.rau), lower=limit, test.ambient=TRUE)
#
e.out

summary(e.out$FDR <= 0.001)


# Concordance by testing with FDR and limited
table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)

sce2 <- sce.rau[,which(e.out$FDR <= 0.001)]


#### Normalize ####

library(scran)
library(scuttle)
library(scater)
clusters <- quickCluster(sce2)
sce2 <- computeSumFactors(sce2, cluster=clusters)


# evaluate ambinat RNA contamination in the empty droplets
amb <- metadata(e.out)$ambient[,1]
head(amb)

sce2 <- logNormCounts(sce2)
#
set.seed(1000)
# modeling variables
dec.pbmc <- modelGeneVarByPoisson(sce2)
# calcualte top features
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)
#
set.seed(1000)
# Evaluate PCs
sce2 <- denoisePCA(sce2, subset.row=top.pbmc, technical=dec.pbmc)
# make UMAP plot
sce2 <- runUMAP(sce2, dimred="PCA")
#
g <- buildSNNGraph(sce2, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce2) <- factor(clust)
#
plotUMAP(sce2,colour_by="label")

#
# remove ambient RNA 
library(scater)
stripped <- sce2[names(amb),]
out <- removeAmbience(counts(stripped), ambient=amb,groups = colLabels(stripped))
#
# load correccted counts into scce object
counts(stripped, withDimnames=FALSE) <- out
stripped <- logNormCounts(stripped)
#
ensmbl_id <- rowData(sce2)$ID[rowData(sce2)$Symbol=="Hba-a1"]
plotExpression(sce2, x="label", colour_by="label", features="Col1a1") +
  ggtitle("Before")

plotExpression(stripped, x="label", colour_by="label", features="Col1a1") + 
  ggtitle("After")

### you left off here
# summary of last things #
# finished doublet detection
# looked at using SoupX, but it needs empty droplets
# looked at https://doi.org/10.1186/s13059-023-02978-x
# looked for other droplet detection methods that work with current data (counts matrix)
# following this tutorial https://rockefelleruniversity.github.io/scRNA-seq/exercises/answers/exercise2_answers.html
# testing on Rau dataset
# will need to run this as a slurm job
# need to decide if empty dropletes/ambient should be run before doublet
# DropletUtils pipeline takes about 15 minutes to run