libs <- c("tidyverse", "Seurat", "SeuratDisk", "Biobase", "MuSiC", 
"reshape2", "SingleCellExperiment") # list libraries here
lapply(libs, function(x){library(x, character.only = T, quietly = T, verbose = F)})
rm(libs)

# Load snRNAseq
sn.anno <- LoadH5Seurat("data/processed/single_cell/celltype_labeled.h5seurat")

# Load phenotype data + markers

phenotypes <- read.csv("data/processed/bulk/pheno_table.csv")
markers <- read.csv("data/processed/single_cell/cluster_markers.csv")

# Load whole bulk RNAseq
bulk.all <- read.csv("data/processed/bulk/all_counts.csv", row.names = 1, check.names = F)

# Run MuSiC with subset genes
bulk.es.exp <- bulk.all[markers$gene,] |> 
  sapply(as.integer) %>%
  ExpressionSet(assayData = .) |> 
  exprs()
row.names(bulk.es.exp) <- markers$gene

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(sn.anno, assay = "RNA")

# Exclude specified clusters
cells <- levels(Idents(sn.anno))

# Subset to measured genes
bulk.es.exp <- bulk.es.exp[!is.na(rowMeans(bulk.es.exp)), ]

# Use MuSiC to estimate cell type proportions
decon <- music_prop(bulk.mtx = bulk.es.exp, sc.sce = sce, markers = markers$gene,
                    clusters = "ident", samples = "orig.ident",
                    select.ct = cells)

# Turn MuSiC output into graph-friendly dataframe
decon.melt <- reshape2::melt(decon$Est.prop.weighted)
colnames(decon.melt) = c('new.id', 'CellType', 'Prop')


# Divide fraction and whole data
decon.frac <- decon.melt|> 
  merge(phenotypes) |> 
  subset(type == "fraction")


decon.whole <- decon.melt|> 
  merge(phenotypes) |> 
  subset(type == "whole")


#Save outputs
if(!dir.exists("data/processed/compositions")){
  dir.create("data/processed/compositions")
}
write.csv(decon.whole, "data/processed/compositions/whole_samples.csv", row.names = F)
write.csv(decon.frac, "data/processed/compositions/fraction_samples.csv", row.names = F)
