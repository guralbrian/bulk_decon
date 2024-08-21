# The aim of this script is to run GO on DEGs that have changed in their signfiicance
# by more than 10x after adjusting for cell type proportions

# ClusterProfiler analysis
libs <- c("tidyverse", "clusterProfiler", "org.Mm.eg.db", "stringr", "data.table") # list libraries here
lapply(libs, require, character.only = TRUE)
rm(libs)

# Get commandArgs
dir.path <- "results/supp_data/DEGs/"
args <- commandArgs(trailingOnly = TRUE)
file_path <- paste0(dir.path, "DEGs_", as.character(args[1]), "_", as.character(args[2]), ".csv")  # Path to the CSV file

# Load the CSV file
data <- read.csv(file_path)

# Function to perform GO enrichment with clusterProfiler
runGo <- function(data, onto){
  # Get a list of significant genes
  sig.genes <- data %>% 
    pull(gene)
  
  # Convert common gene names to ENSEMBL IDs for clusterProfiler
  gene.df <- bitr(sig.genes, fromType = "SYMBOL",
                  toType = c("ENSEMBL"),
                  OrgDb = org.Mm.eg.db)
  
  # Check for enrichment within biological process gene clusters
  ego <- enrichGO(gene          = gene.df$ENSEMBL,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = onto,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.1,
                  qvalueCutoff  = 1,
                  readable      = TRUE) 
  return(ego)
}

# Perform GO enrichment analysis
go.df <- runGo(data, "BP")

# make into df
go.df <- go.df@result |> 
           mutate(qscore = -log(p.adjust, base=10))

# Make data directory to store CSVs
data.path <- "results/supp_data/go_terms/"
if(!dir.exists(data.path)){
  dir.create(data.path)
}

# Save the result
output_file <- paste0(data.path, basename(file_path))
write.csv(go.df, output_file)