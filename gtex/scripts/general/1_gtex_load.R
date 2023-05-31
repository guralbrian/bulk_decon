library(dplyr)
libs <- c("Seurat", "zellkonverter", "curl", "SeuratDisk")
lapply(libs, require, character.only = T)
setwd("/proj/raulab/users/brian/r_projects/gtex")

#### initial data processing and organization #####

# GTEx snRNAseq
url <- "https://storage.googleapis.com/gtex_analysis_v9/snrna_seq_data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad" 
temp <- tempfile(fileext = ".h5ad")
curl::curl_download(url, temp)
gtex.sn <- zellkonverter::readH5AD(temp, verbose=F, layers=T, varm=F, obsm=F, varp=F, obsp=F, uns=F)
gtex.sn <- as.Seurat(
  gtex.sn,
  counts = "counts",
  data = "X",
  assay = NULL,
  project = "SingleCellExperiment"
)

# clean up metadata
colnames(gtex.sn@meta.data)[c(2,3,11)] <- c("nCount_RNA", "nFeature_RNA", "PercentMito")
gtex.sn <- RenameAssays(gtex.sn, "originalexp" = "RNA")  |>
  subset(subset = tissue == "heart")

meta <- names(gtex.sn@meta.data)[sapply(gtex.sn@meta.data, is.factor)]

for(i in meta){
 gtex.sn@meta.data[[i]] <- droplevels(gtex.sn@meta.data[[i]])
}

SaveH5Seurat(gtex.sn, "data/processed/internal/sn_gtex_lv_match.h5seurat", overwrite = T)

# GTEx bulk
url <- "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_heart_left_ventricle.gct.gz"
temp <- tempfile(fileext = ".gct.gz")
curl::curl_download(url, temp)
gtex.bk <- read.delim(temp, 
                      skip = 2, row.names = 1, stringsAsFactors = F)

# sum expression of gene isoforms
gtex.bk <- gtex.bk[,-1] |>
  group_by(Description) |>
  summarise(across(everything(), sum)) |>
  as.data.frame()

# rename columns to match Participant.ID in gtex.sn
names <- colnames(gtex.bk)[grepl("GTEX", colnames(gtex.bk))] |>
  strsplit(colnames(gtex.bk), split = "[.]") |>
  lapply('[', 2) 

names <- paste0("GTEX-", names)

colnames(gtex.bk)  <- c("gene", names)

write.csv(gtex.bk, "data/processed/internal/gtex_lv_counts_summed.csv")
