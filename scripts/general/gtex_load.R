library(dplyr)
libs <- c("Seurat", "zellkonverter", "curl")
require(libs)
lapply(libs, require, character.only = T)

#### initial data processing and organization #####

# GTEx snRNAseq
url <- "https://storage.googleapis.com/gtex_analysis_v9/snrna_seq_data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad" 
temp <- tempfile(fileext = ".h5ad")
curl::curl_download(url, temp)
gtex.sn <- zellkonverter::readH5AD(temp, verbose=T, layers=T, varm=F, obsm=F, varp=F, obsp=F, uns=F)
gtex.sn <- as.Seurat(
  gtex.sn,
  counts = "counts",
  data = "X",
  assay = NULL,
  project = "SingleCellExperiment"
)
colnames(gtex.sn@meta.data)[c(2,3,11)] <- c("nCount_RNA", "nFeature_RNA", "PercentMito")
gtex.sn <- RenameAssays(gtex.sn, "originalexp" = "RNA")  |>
              subset(subset = tissue == "heart")
gtex.sn$Participant.ID <- droplevels(gtex.sn$Participant.ID)

SaveH5Seurat(gtex.sn, "data/processed/internal/sn_gtex_lv_match.h5seurat", overwrite = T)


# GTEx bulk
url <- "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_heart_left_ventricle.gct.gz"
temp <- tempfile(fileext = ".gct.gz")
curl::curl_download(url, temp)

bk_gtex <- read.delim(temp, 
  skip = 2, row.names = 1, stringsAsFactors = F)

bk_gtex <- bk_gtex[,-1] |>
  group_by(Description) |>
  summarise(across(everything(), sum)) |>
  as.data.frame()

write.csv(bk_gtex, "data/processed/internal/gtex_lv_counts_summed.csv", overwrite = F)
