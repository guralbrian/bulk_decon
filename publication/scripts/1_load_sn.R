libs <- c("Seurat", "zellkonverter", "curl", "SeuratDisk",
          "TabulaMurisSenisData", "tidyverse", "Matrix")
lapply(libs, require, character.only = T, lib.loc = "tools/r_libs/common")


# Rau internal (raw)####

rau.raw.dir <- "data/raw/rau_sn/patterson/raw_matrices"
rau.samples <- list.files(rau.raw.dir)

# List the .mtx file and read
LoadRau <- function(sample){
  rau.samp.dir <- paste0(rau.raw.dir,"/", sample)
  rau.files <- list.files(rau.samp.dir)
  mtx.rau <- rau.files[[grep("mtx", rau.files)]]
  mtx.rau <- readMM(paste0(rau.samp.dir, "/",mtx.rau))
  
  # Get list the barcodes .tsv file
  
  bar.rau <- rau.files[[grep("barcodes", rau.files)]]
  bar.rau <- readLines(paste0(rau.samp.dir, "/",bar.rau))
  
  # Get list the genes .tsv file
  ft.rau <- rau.files[[grep("features", rau.files)]]
  ft.rau <- read.table(paste0(rau.samp.dir, "/",ft.rau))
  
  colnames(mtx.rau) <- paste0(bar.rau, "-rau")
  rownames(mtx.rau) <- ft.rau$V2
  
  sn.rau <- CreateSeuratObject(mtx.rau)
  sn.rau$origin <- "brian"
  sn.rau$orig.ident <- paste0("brian_", sample)
  sn.rau$PercentMito <- PercentageFeatureSet(sn.rau, pattern = "^mt-")
  sn.rau$PercentRibo <- PercentageFeatureSet(sn.rau, pattern = "^Rpl|^Rps")
  return(sn.rau)
}

rau.list <- lapply(rau.samples, function(x){LoadRau(x)})

SaveH5Seurat(rau.list[[1]], "data/raw/rau_sn/patterson/seurat/b6_1")
SaveH5Seurat(rau.list[[2]], "data/raw/rau_sn/patterson/seurat/b6_2")
