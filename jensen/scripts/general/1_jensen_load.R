libs <- c("Seurat", "zellkonverter", "curl", "SeuratDisk",
          "TabulaMurisSenisData", "tidyverse")
lapply(libs, require, character.only = T)

#### initial data processing and organization #####
# Tabula Muris ####
sn_muris <- TabulaMurisSenisDroplet(
  tissues = "Heart_and_Aorta",
  processedCounts = FALSE,
  reducedDims = TRUE,
  infoOnly = FALSE)[[1]] 

sn_muris_seurat <- sn_muris |>
  SummarizedExperiment::assay("counts")|>
  as.matrix() |>
  Seurat::CreateSeuratObject()

# Add metadata from SingleCellExperiment to Seurat
for(i in colnames(colData(sn_muris))){
  sn_muris_seurat <- Seurat::AddMetaData(object = sn_muris_seurat, col.name = i, metadata = colData(sn_muris)[[i]])
}

# Tabula Muris didn't map to mitochondrial genes
# but we need the data slot for later
#! To publish, I'd need to remap their data
sn_muris_seurat$PercentMito <- 0
sn_muris_seurat$origin <- "tabula_muris"

SaveH5Seurat(sn_muris_seurat, "jensen/data/processed/tabula_muris")

# Wu 2021 ####

# Import data from Wu 2021 doi.org/10.1152/physiolgenomics.00016.2021
geo <- "GSM5471468"
filePaths = getGEOSuppFiles("GSM5471468")

# Unzip files
paths <- rownames(filePaths)
for(i in paths){
  gunzip(i)
}

# Read and format as Seurat
wu.mt <- readMM("GSM5471468/GSM5471468_matrix.mtx")
cell.ids <- readLines("GSM5471468/GSM5471468_barcodes.tsv")
genes <- read.table("GSM5471468/GSM5471468_genes.tsv")
colnames(wu.mt) <- paste0(cell.ids, "-wu")
rownames(wu.mt) <- genes$V2
sn.wu <- CreateSeuratObject(wu.mt)
sn.wu$origin <- "wu"
sn.wu$orig.ident <- "wu"
SaveH5Seurat(sn.wu, "jensen/data/processed/single_cell/wu_2021")


# Martini 2019 ####

# Import data from Martini 2019 doi.org/10.1161/CIRCULATIONAHA.119.041694
geo <- "GSE122930"
dir <- "./jensen/data/processed/geo"
filePaths = getGEOSuppFiles(geo, baseDir = dir)
paths <- rownames(filePaths)

# Unzip files
for(i in paths){
  if(endsWith(i, '.gz')){
    gunzip(i)
  }
}

# select files ending in .mtx or .tsv
files <- list.files(paste0(dir, "/", geo))
files.mtx <- files[endsWith(files, '.mtx')]
files.tsv <- files[endsWith(files, '.tsv')]

# group by lapply(names, strsplit( "_"), "[[", 2)

treatment <- lapply(strsplit(files.mtx, "_"), "[[", 2)
timepoint <- lapply(strsplit(files.mtx, "_"), "[[", 3)
samples.martini <- paste0(treatment,"_", timepoint)

# read in each file and merge

for(i in samples.martini){
  # List the .mtx file and read
  mtx.file <- files.mtx[grepl(i, files.mtx)]
  mtx <- readMM(paste0(dir, "/",geo, "/", mtx.file))
  
  # Get list the barcodes .tsv file
  tsv.bar <- files.tsv[grepl(i, files.tsv) & grepl("barcodes", files.tsv)]
  tsv.bar <- readLines(paste0(dir, "/",geo, "/", tsv.bar))
  
  # Get list the genes .tsv file
  tsv.genes <- files.tsv[grepl(i, files.tsv) & grepl("genes", files.tsv)]
  tsv.genes <- read.table(paste0(dir, "/",geo, "/", tsv.genes))
  
  # Format and return sparse matrix
  colnames(mtx) <- paste0(tsv.bar, "-mar")
  rownames(mtx) <- tsv.genes$V2
  assign(i, mtx)
}

# Convert to Seurats
for(i in 1:length(samples.martini)){
  seurat <- CreateSeuratObject(get(samples.martini[i]))
  seurat$martini.cond <- samples.martini[i]
  assign(samples.martini[i], seurat)
}

# Merge Seurats and save
sn.martini <- merge(get(samples.martini[1]), get(samples.martini[2])) |>
  merge(get(samples.martini[3])) |>
  merge(get(samples.martini[4]))
sn.martini$origin <- "martini"
sn.martini$orig.ident <- paste0("martini_", sn.martini$martini.cond)
SaveH5Seurat(sn.martini, "jensen/data/processed/single_cell/martini_2019")

# Froese ####

# Function specific to loading data from links in Froese paper
GetFroese <-function(url){
  GET(url, write_disk(tf <- tempfile(fileext = ".xlxs")))
  df <- read_excel(tf)
  # make just gene, log2cpm (by sham sample)
  froese <- df[,c(3,14:16)] |>
    as.data.frame()
  froese <- froese[!duplicated(froese$SYMBOL),]
}

# load each fraction
froese_cm <- GetFroese("https://ars.els-cdn.com/content/image/1-s2.0-S2589004222002358-mmc2.xlsx")
froese_fb <- GetFroese("https://ars.els-cdn.com/content/image/1-s2.0-S2589004222002358-mmc4.xlsx")
froese_ec <- GetFroese("https://ars.els-cdn.com/content/image/1-s2.0-S2589004222002358-mmc6.xlsx")

# Merge fractions
froese <- Reduce(function(x,y) merge(x,y,by="SYMBOL",all=TRUE), 
                 list(froese_cm,froese_fb,froese_ec))
rownames(froese) <- froese$SYMBOL
froese <- froese[,-1]
rm(froese_cm,froese_fb,froese_ec)
write.csv(froese, "jensen/data/processed/froese_fractions")

# Christoph fractions + Jensen whole df + Froese fractions, in CPM
jensen_whole <- read.csv("jensen/data/raw/jensen_counts.csv", header = 1, row.names = 1)
rau_fractions <- read.csv("jensen/data/processed/celltype_counts.csv", row.names = 1)

# Merge non-cpm data
all_bulk <- merge(rau_fractions, jensen_whole, by = "row.names")
rownames(all_bulk) <- all_bulk$Row.names
all_bulk <- all_bulk[,-1]

# Calculate the library size for each sample
lib_sizes <- matrixStats::colSums2(as.matrix(all_bulk))

# Normalize counts by library size and multiply by 1,000,000
cpm <- t(apply(all_bulk, 1, function(x) (x / lib_sizes) * 1e6)) 
cpm <- cpm + 0.0000001
log2_cpm <- log2(cpm)

# Merge all data
all_bulk <- merge(log2_cpm, froese, by = "row.names")
rownames(all_bulk) <- all_bulk$Row.names
all_bulk <- all_bulk[,-1]

write.csv(all_bulk, "jensen/data/processed/jensen_rau_froese_cpm")

# Phenotype/metadata for bulk RNAseq ####

pheno <- read.csv("jensen/data/raw/bulk_phenotypes.csv", header = 1, row.names = 1)
pheno$genotype <- "wt"
pheno[grepl("A", pheno$id), 3] <- "ko"
pheno$treatment <- "control"
pheno[grepl("Lx", pheno$id), 4] <- "mi"

all_pheno <- data.frame(id = c(colnames(all_bulk)),
                        type = c(rep("fraction_christoph", length(fractions_pheno$colname)),
                                 rep("whole_jensen", length(rownames(pheno))),
                                 rep("fraction_froese", length(colnames(froese)))
                        )
)

rownames(all_pheno) <- c(rownames(pheno), fractions_pheno$colname, colnames(froese))

# Name cell types by origin-specific nomenclature
all_pheno$fraction <- "whole_tissue"
all_pheno[grepl("CM", all_pheno$id), 3] <- "Cardiomyocytes"
all_pheno[grepl("Fib", all_pheno$id), 3] <- "Fibroblasts"
all_pheno[grepl("FB", all_pheno$id), 3] <- "Fibroblasts"
all_pheno[grepl("Endo", all_pheno$id), 3] <- "Endothelial Cells"
all_pheno[grepl("EC_", all_pheno$id), 3] <- "Endothelial Cells"

all_pheno$origin <- "Jensen"
all_pheno[grepl("B6", all_pheno$id), 4] <- "Rau"
all_pheno[grepl("Sham", all_pheno$id), 4] <- "Froese"

write.csv(all_pheno, "jensen/data/processed/jensen_rau_froese_pheno")