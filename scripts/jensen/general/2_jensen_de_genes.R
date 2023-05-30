# Pseudobulk-bulk DE filtering on GTEx samples
# Compare subject matched bulk RNAseq to pseudobulk single-nucleus (sn) RNAseq

# load libraries 
library(dplyr) # load first since we'll need dplyr::lapply()
libs <- c("Seurat", "DESeq2", "SeuratDisk") # list libraries here
lapply(libs, require, character.only = T)

source("scripts/functions/decon_all.R") # load our own functions

# load data
gtex.sn <- LoadH5Seurat("data/processed/internal/sn_gtex_lv_match.h5seurat")
gtex.bk <- read.csv("data/processed/internal/gtex_lv_counts_summed.csv", check.names = F, row.names = 1)

# keep only overlapping subjects
gtex.bk.sub <- gtex.bk[,levels(gtex.sn$Participant.ID)]
rownames(gtex.bk.sub) <- gtex.bk$gene

# find DE genes
good.genes <- FilterBulkSingleNucleus(gtex.sn, 
                                     gtex.bk.sub,
                                     min.rna.ft = 200,
                                     max.rna.ft = 3000,
                                     min.rna.ct = 800,
                                     max.mt.pt = 0.01,
                                     max.rb.pt = 0.01,
                                     scrublet.score = 0.4,
                                     group = "Participant.ID",
                                     lfc.thresh = log2(5),
                                     change = "greaterAbs") 

# Make new seurat with stable genes 
counts <- GetAssayData(gtex.sn, slot="counts", assay="RNA")[good.genes,] |>
  CreateSeuratObject()

# add meta data
meta.features <- colnames(gtex.sn@meta.data)
for(i in meta.features){
  counts <- AddMetaData(counts, gtex.sn@meta.data[[i]], col.name = i)
}

# reassign and clear junk
gtex.sn <- counts
rm(counts)
gc()

SaveH5Seurat(gtex.sn, "data/processed/internal/sn_gtex_no_de.h5seurat", overwrite = T)