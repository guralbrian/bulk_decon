https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5471468&format=file&file=GSM5471468%5Fmatrix%2Emtx%2Egz
library(GEOquery)
library(Matrix)
geo <- "GSM5471468"
filePaths = getGEOSuppFiles("GSM5471468")

paths <- rownames(filePaths)
for(i in paths){
  gunzip(i)
}

library(Matrix)
files <- list.files()
wu.mt <- readMM("GSM5471468/GSM5471468_matrix.mtx")

colName(wu.mt)
gsub(genes[[1]], pattern = "t", "__")
strsplit(genes, "t")

cell.ids <- readLines("GSM5471468/GSM5471468_barcodes.tsv")
genes <- read.table("GSM5471468/GSM5471468_genes.tsv")
colnames(wu.mt) <- paste0(cell.ids, "-wu")
rownames(wu.mt) <- genes$V2

length(unique(genes$V2))
length(genes$V2)
sn.wu <- CreateSeuratObject(wu.mt)
sn.wu$origin <- "wu"
sn.wu$orig.ident <- "wu"

sn <- LoadH5Seurat("jensen/data/processed/single_cell/05102023.h5seurat")

sn.all <- merge(sn, sn.wu)
sn.all$PercentMito <- PercentageFeatureSet(sn.all, pattern = "^mt-")
sn.all$origin <- as.factor(sn.all$origin)
unique(sn.all$orig.ident)
sn.clust.1 <- ClusterSeurat(sn.all, max.rna.ft = 10000,
                            max.mt.pt = 1, res = 0.05, regress.by = c("orig.ident"))

dim.ft  <- c("seurat_clusters", "origin", "orig.ident")
feat.ft <- c("PercentMito", "nCount_RNA")
vln.ft <- c("PercentMito")
vln.groups <- c("seurat_clusters")

plotUMAP(data = sn.clust.1,
         dim.ft = dim.ft,
         feat.ft = feat.ft,
         nrow = 3, 
         ncol = 2)

do_DimPlot(sn.clust.1, group.by = "orig.ident")
