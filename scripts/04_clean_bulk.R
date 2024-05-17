libs <- c("tidyverse",  "reshape2",  "makeunique") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#Load bulk RNAseq dataset

all_bulk <-  read.csv("data/processed/bulk/all_bulk_gene.csv", row.names = 1)

phenotypes <- data.frame(id = colnames(all_bulk))

# Add origin column
phenotypes <- phenotypes |> 
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
      str_detect(id, "Lx") ~ "MI",
      str_detect(id, "B6_") ~ NA,
      .default = "Sham")
  )
phenotypes[which(phenotypes$id == "Ako_Lx4"), "genotype"] <-"WT"
phenotypes[which(phenotypes$id == "wt_Lx3"), "genotype"] <-"cmAKO"

# Reformat names 

phenotypes <- phenotypes |> 
  mutate(
    gene_treat = paste(genotype, treatment)
  )
phenotypes$gene_treat <- make_unique(phenotypes$gene_treat)

phenotypes <- phenotypes |>  
  mutate(
    new.id = case_when(
      type == "whole" ~ gene_treat,
      type == "fraction" ~ id
    )
  )

colnames(all_bulk) <- phenotypes$new.id

# Save
write.csv(all_bulk, "data/processed/bulk/all_counts.csv", row.names = T)
write.csv(phenotypes, "data/processed/bulk/pheno_table.csv", row.names = F)


