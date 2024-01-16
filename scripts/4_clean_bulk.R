libs <- c("tidyverse",  "reshape2",  "makeunique") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)


#Load cell type fraction bulk RNAseq dataset
fractions <- read.csv("data/raw/rau_fractions/celltype_counts.csv", row.names = 1)
fractions.pheno <- read.csv("data/raw/rau_fractions/celltype_pheno.csv")

# Load whole bulk RNAseq
jensen_bulk <- readxl::read_xlsx("data/raw/jensen/jensen_counts_correct.xlsx")[,-c(20,21)]|>
  as.data.frame()

# Reformat names and such
# Will be changed after in-house alignment pipeline is complete
colnames(jensen_bulk)[c(20,21)] <- c("WT MI", "KO MI")

phenotypes_real <- jensen_bulk[1,] |> 
  pivot_longer(cols = -c(1:5), 
               names_to = "gene_treat", 
               values_to = "id") |>
  mutate(gene_treat = str_remove(gene_treat, "\\..*")) |>
  select(gene_treat, id)

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
write.csv(bulk.all, "data/processed/bulk/all_counts.csv", row.names = F)
write.csv(phenotypes_real, "data/processed/bulk/jensen_pheno.csv", row.names = F)
write.csv(jensen_bulk, "data/processed/bulk/jensen_bulk_clean.csv")
