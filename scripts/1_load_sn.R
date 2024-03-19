# This script is meant to merge barcodes, features, and counts into a Seurat object
# The data is from Drs. Mikayla Patterson and Christoph Rau. 
# It is snRNAseq from male C57B6/J left ventricles. There are two datasets, each
# prepared from the composite of 4-5 mouse isolations

# Load Libraries
libs <- c("Seurat", "Matrix", "SeuratDisk")
lapply(libs, require, character.only = T)

# Get commandArgs
args <- commandArgs(trailingOnly = TRUE)
print(args)
sample_name <- as.character(args[1])
print(sample_name)

path <- paste0("data/raw/single_cell/", sample_name) 

rau.files <- path |> list.files()

# read .mtx, barcodes, and features files + merge into Seurat

  # Get the counts matrix
  mtx.rau <- rau.files[[grep("mtx", rau.files)]]
  mtx.rau <- Matrix::readMM(paste0(path, "/",mtx.rau))
  print("loaded matrix")
  # Get list the barcodes .tsv file
  bar.rau <- rau.files[[grep("barcodes", rau.files)]]
  bar.rau <- readLines(paste0(path, "/",bar.rau))
  print("loaded barcodes")  
  # Get list the genes .tsv file
  ft.rau <- rau.files[[grep("features", rau.files)]]
  ft.rau <- read.table(paste0(path, "/",ft.rau))
  print("loaded genes")  
  # Format to show origin for later if needed
  colnames(mtx.rau) <- paste0(bar.rau, "-rau")
  rownames(mtx.rau) <- ft.rau$V2 |> make.unique()
  
  sn.rau <- Seurat::CreateSeuratObject(mtx.rau, min.cells = 3, min.features = 10)
  sn.rau$origin <- "brian"
  sn.rau$orig.ident <- paste0("brian_", sample_name)
  print("made seurat")
  # Measure mitochondrial and ribosomal gene percentages 
  sn.rau$PercentMito <- Seurat::PercentageFeatureSet(sn.rau, pattern = "^mt-")
  sn.rau$PercentRibo <- Seurat::PercentageFeatureSet(sn.rau, pattern = "^Rpl|^Rps")

# Save
SeuratDisk::SaveH5Seurat(sn.rau, paste0("data/processed/single_cell/unprocessed/",sample_name),
                         overwrite = T,
                         verbose = T)
