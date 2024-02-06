# This script is meant to determine if any of the nuclei in our snRNAseq dataset are activated fibroblasts

# Load libs and data
BiocManager::install("gt")
libs <- c("tidyverse", "Seurat", "SeuratDisk",  "Biobase", "gt") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)


sn <- LoadH5Seurat("data/processed/single_cell/merged_no_doublets.h5seurat")

# Exclude CMs + get ratios
sn.mark$cell_type <- Idents(sn.mark) 
cell_counts <- sn.mark@meta.data |> 
  subset(cell_type != "Cardiomyocytes") |> 
  pull("cell_type") |>
  table()
cell_ratios <- cell_counts / sum(cell_counts) 

# Format for tables
cell_ratios <- cell_ratios |> 
  as.data.frame()
colnames(cell_ratios) <- c("Cell type", "Proportion")
cell_ratios |> 
  mutate(
    Proportion = Proportion * 10^3,
    Proportion = floor(Proportion)
  ) |> 
  filter(
    `Cell type` != "Cardiomyocytes") |> 
  gt() |>
  tab_header(
    title = "Nuclei per 1000"
  )
