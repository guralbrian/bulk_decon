# This script is meant to merge barcodes, features, and counts into a Seurat object
# The data is from Drs. Mikayla Patterson and Christoph Rau. 
# It is snRNAseq from male C57B6/J left ventricles. There are two datasets, each
# prepared from the composite of 4-5 mouse isolations

# Load Libraries
libs <- c("Seurat", "Matrix", "SeuratDisk")
lapply(libs, require, character.only = T)

rau.raw.dir <- "data/raw/single_cell"
rau.samples <- list.files(rau.raw.dir)

# read .mtx, barcodes, and features files + merge into Seurat
LoadRau <- function(sample){
  # Specify the sample directory
  rau.samp.dir <- paste0(rau.raw.dir,"/", sample)
  rau.files <- list.files(rau.samp.dir)
  
  # Get the counts matrix
  mtx.rau <- rau.files[[grep("mtx", rau.files)]]
  mtx.rau <- Matrix::readMM(paste0(rau.samp.dir, "/",mtx.rau))
  
  # Get list the barcodes .tsv file
  bar.rau <- rau.files[[grep("barcodes", rau.files)]]
  bar.rau <- readLines(paste0(rau.samp.dir, "/",bar.rau))
  
  # Get list the genes .tsv file
  ft.rau <- rau.files[[grep("features", rau.files)]]
  ft.rau <- read.table(paste0(rau.samp.dir, "/",ft.rau))
  
  # Format to show origin for later if needed
  colnames(mtx.rau) <- paste0(bar.rau, "-rau")
  rownames(mtx.rau) <- ft.rau$V2
  
  sn.rau <- Seurat::CreateSeuratObject(mtx.rau)
  sn.rau$origin <- "brian"
  sn.rau$orig.ident <- paste0("brian_", sample)
  
  # Measure mitochondrial and ribosomal gene percentages 
  sn.rau$PercentMito <- Seurat::PercentageFeatureSet(sn.rau, pattern = "^mt-")
  sn.rau$PercentRibo <- Seurat::PercentageFeatureSet(sn.rau, pattern = "^Rpl|^Rps")
  return(sn.rau)
}

# Apply to each sample directory
rau.list <- lapply(rau.samples, function(x){LoadRau(x)})

# Save
SeuratDisk::SaveH5Seurat(rau.list[[1]], "data/processed/single_cell/b6_1")
SeuratDisk::SaveH5Seurat(rau.list[[2]], "data/processed/single_cell/b6_2")
