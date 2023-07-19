library(dplyr)
libs <- c("Seurat", "zellkonverter", "curl", "SeuratDisk", "tidyverse", "Matrix", "data.table")
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

SaveH5Seurat(gtex.sn, "gtex/data/processed/internal/sn_gtex_lv_match.h5seurat", overwrite = T)

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

names <- gsub("[.]", "-", colnames(gtex.bk))

colnames(gtex.bk)  <- c("gene", names[-1])

write.csv(gtex.bk, "gtex/data/processed/internal/gtex_lv_counts_summed.csv")



#### Extra sn datasets

## Simonson 2023
# get list of files
sim.dir <- "gtex/data/raw/external/single_cell/simonson_2023/"
sim.files <- list.files(sim.dir)

# load in files
sim.cells <- readLines(paste0(sim.dir,sim.files[[1]]))
sim.genes <- read.table(paste0(sim.dir,sim.files[[2]]))
sim.counts <- readMM(paste0(sim.dir,sim.files[[3]]))
sim.meta <-  fread(paste0(sim.dir,sim.files[[4]]), header = T)

# format meta data classes
data_types <- slice(sim.meta, 1)

# Remove the first row from the data
sim.meta <- slice(sim.meta, -1)

# Get the column names for each type
group_cols <- names(data_types)[data_types == "group"]
numeric_cols <- names(data_types)[data_types == "numeric"]

# Convert columns to their appropriate data types
sim.meta <- sim.meta %>%
  mutate(across(all_of(group_cols), as.factor))    |>
  mutate(across(all_of(numeric_cols), as.numeric)) |>
  mutate(NAME = paste0(NAME, "-sim"))

rownames(sim.meta) <- sim.meta$NAME
sim.meta <- sim.meta[,-1]

colnames(sim.meta)[which(colnames(sim.meta) == "donor_id")] <- "Participant.ID"
colnames(sim.meta)[which(colnames(sim.meta) == "n_genes")] <- "nFeature_RNA"
colnames(sim.meta)[which(colnames(sim.meta) == "n_umi")] <- "nCount_RNA"
colnames(sim.meta)[which(colnames(sim.meta) == "cellranger_percent_mito")] <- "PercentMito"


## Merge into seurat

colnames(sim.counts) <- paste0(sim.cells, "-sim")
rownames(sim.counts) <- sim.genes$V2
sim.sn <- CreateSeuratObject(sim.counts)


# Add meta data
for(i in colnames(sim.meta)){
  sim.sn <- Seurat::AddMetaData(object = sim.sn, col.name = i, metadata = sim.meta[[i]])
}

SaveH5Seurat(sim.sn, paste0(sim.dir, "simonson2023_raw"))



#### Purified Cell Types #### 
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE218251&format=file&file=GSE218251%5FProcessedData%2Exlsx

# LV Fibroblasts
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218251

# Heart ECs
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114607

# Bulk on nuclei of CMs
https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.120.051921#d1e1048

