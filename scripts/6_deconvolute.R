
libs <- c("tidyverse", "Seurat", "SeuratDisk", "SCpubr", "Biobase", "MuSiC", 
"reshape2", "SingleCellExperiment") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load snRNAseq
sn <- LoadH5Seurat("data/processed/single_cell/celltype_labeled.h5seurat")

# Load phenotype data + markers
fractions.pheno <- read.csv("data/raw/rau_fractions/celltype_pheno.csv")
phenotypes_real <- read.csv("data/processed/bulk/jensen_pheno.csv")
all.markers <- read.csv("data/processed/single_cell/cluster_markers.csv")

# Load whole bulk RNAseq
bulk.all <- read.csv("data/processed/bulk/all_counts.csv", row.names = 1, check.names = F)

# Run MuSiC with subset genes
bulk.es.exp <- bulk.all[all.markers$gene,] |> 
  sapply(as.integer) %>%
  ExpressionSet(assayData = .) |> 
  exprs()
row.names(bulk.es.exp) <- all.markers$gene


# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(sn, assay = "RNA")

# Exclude specified clusters
cells <- levels(Idents(sn))

# Use MuSiC to estimate cell type proportions
decon <- music_prop(bulk.mtx = bulk.es.exp, sc.sce = sce, markers = all.markers$gene,
                    clusters = "ident", samples = "orig.ident",
                    select.ct = cells)

# Turn MuSiC output into graph-friendly dataframe
decon.melt <- reshape2::melt(decon$Est.prop.weighted)
colnames(decon.melt) = c('Sub', 'CellType', 'Prop')

# Divide fraction and whole data
#decon.melt$version <- "music" 
#decon.melt <- rbind(decon.melt, props.bisque)
colnames(fractions.pheno)[4] <- "Sub"
decon.frac <- decon.melt|> 
  subset(Sub %in% fractions.pheno$Sub) |> 
  merge(fractions.pheno)


decon.whole <- decon.melt|> 
  #mutate(Sub = paste0("S", Sub)) |> 
  subset(Sub %in% phenotypes_real$de_id) |> 
  mutate(Genotype = factor(case_when(
    str_detect(Sub, "WT") ~ "WT",
    str_detect(Sub, "KO") ~ "KO")),
    Treatment = factor(case_when(
      str_detect(Sub, "sham") ~ "Sham",
      str_detect(Sub, "MI") ~ "TAC")))


#Save outputs
write.csv(decon.whole, "data/processed/compositions/whole_samples.csv", row.names = F)
write.csv(decon.frac, "data/processed/compositions/fraction_samples.csv", row.names = F)
