libs <- c("tidyverse",  "reshape2",  "makeunique") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#Load bulk RNAseq dataset

#fractions.pheno <- read.csv("data/raw/rau_fractions/celltype_pheno.csv")
all_bulk <-  read.csv("data/processed/bulk/all_bulk_gene.csv", row.names = 1)

# Load whole bulk RNAseq
#jensen_bulk <- readxl::read_xlsx("data/raw/jensen/jensen_counts_correct.xlsx")[,-c(20,21)]|>
#  as.data.frame()

phenotypes <- data.frame(id = colnames(all_bulk))

# Add origin column
phenotypes |> 
  mutate(
    type = case_when(
      str_detect(id, "B6_") ~ "fraction",
      .default = "whole"),
    cell.type = case_when(
      str_detect(id, "CM_") ~ "Cardiomyocyte",
      str_detect(id, "Fib_") ~ "Fibroblast",
      str_detect(id, "Endo") ~ "Endothelial Cell",
      .default = NA),
    sex = case_when(
      str_detect(id, "_M_") ~ "Male",
      str_detect(id, "_F_") ~ "Female",
      .default = NA),
    sex = case_when(
      str_detect(id, "_M_") ~ "Male",
      str_detect(id, "_F_") ~ "Female",
      .default = NA),
    genotype = case_when(
      str_detect(id, "WT") | str_detect(id, "Ako_") ~ "WT",
      str_detect(id, "AKO") | str_detect(id, "wt_")~ "cmAKO",
      .default = NA),
    treatment = case_when(
      str_detect(id, "Lx") ~ "CAD",
      str_detect(id, "B6_") ~ NA,
      .default = "Sham")
  )
# Reformat names and such
# Will be changed after in-house alignment pipeline is complete
colnames(jensen_bulk)[c(20,21)] <- c("WT MI", "KO MI")

phenotypes_real <- jensen_bulk[1,] |> 
  pivot_longer(cols = -c(1:5), 
               names_to = "gene_treat", 
               values_to = "id") |>
  mutate(gene_treat = str_remove(gene_treat, "\\..*")) |>
  dplyr::select(gene_treat, id)

jensen_bulk <- jensen_bulk[-1,-c(1:4)]
colnames(jensen_bulk) <- c("gene", make_unique(unlist(phenotypes_real$gene_treat)))

phenotypes_real$de_id <- colnames(jensen_bulk)[-1]

# Add a suffix to duplicate gene names to make them unique
make.names.jensen_bulk <- function(jensen_bulk){
  duplicated.names <- duplicated(jensen_bulk$gene)
  suffix <- cumsum(duplicated.names)
  jensen_bulk$gene[duplicated.names] <- paste0(jensen_bulk$gene[duplicated.names], ".", suffix[duplicated.names])
  return(jensen_bulk)
}

jensen_bulk <- make.names.jensen_bulk(jensen_bulk)

# Now set the row names
row.names(jensen_bulk) <- jensen_bulk$gene

jensen_bulk <- jensen_bulk[,-1]

jensen_bulk <- mutate_all(jensen_bulk, function(x) round(as.numeric(as.character(x)), digits = 0)) # round to integers
phenotypes_real$de_id <- colnames(jensen_bulk)


common.genes <- row.names(fractions)[row.names(fractions) %in% row.names(jensen_bulk)]

bulk.all <- cbind(jensen_bulk[common.genes,], fractions[common.genes,])

# Save datasets
write.csv(bulk.all, "data/processed/bulk/all_counts.csv", row.names = T)
write.csv(phenotypes_real, "data/processed/bulk/jensen_pheno.csv", row.names = F)
write.csv(jensen_bulk, "data/processed/bulk/jensen_bulk_clean.csv")


# Make fractions phenotypes file
# list files
files <- list.files("data/raw/fastq")[str_detect(list.files("data/raw/fastq"), "B6_")]

# Split name details
files.split <- str_split(files, "_")

# Compile in df
files.df <- data.frame(Sub = files)
files.df$sex <- lapply(files.split, "[[", 2) |> unlist()
files.df$cell.type <- lapply(files.split, "[[", 3) |> unlist()
files.df$replicate <- lapply(files.split, "[[", 5) |> unlist()

# Save
write.csv(files.df, "data/raw/rau_fractions/celltype_pheno.csv", row.names = F)
