# ClusterProfiler analysis
libs <- c("tidyverse", "clusterProfiler", "DESeq2", 
          "org.Mm.eg.db", "viridis", "stringr", "patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Get commandArgs
args <- commandArgs(trailingOnly = TRUE)
model.type <-  as.character(args[1])
#model.type <- "adjusted"
# Load the results and expression matrices
res <- readRDS(paste0("data/processed/models/", model.type,"_de_interaction.RDS")) 
names <- resultsNames(res)[-1]

# make a list of each result
res <- lapply(names, function(x){
  results(res, name=x) |> 
    as.data.frame()
})

names(res) <- names

# Function to perform GO enrichment with clusterProfiler
runGo <- function(data, onto){
# Get a list of significant genes
sig.genes <- data |> 
  filter(padj < 0.1 & abs(log2FoldChange) > 0.583) |> 
  row.names()

# Convert common gene names to ENSEMBLE IDs for clusterProfiler
gene.df <- bitr(sig.genes, fromType = "SYMBOL",
                toType = c("ENSEMBL"),
                OrgDb = org.Mm.eg.db)

gene.list <- bitr(row.names(data), fromType = "SYMBOL",
                toType = c("ENSEMBL"),
                OrgDb = org.Mm.eg.db)

# Check for enrichment within biological process gene clusters
ego <- enrichGO(gene          = gene.df$ENSEMBL,
                OrgDb         = org.Mm.eg.db,
                universe      = gene.list$ENSEMBL,
                keyType       = 'ENSEMBL',
                ont           = onto,
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE) #|> 
        #clusterProfiler::simplify(cutoff = 0.6)
}

# Apply the function
go <- lapply(res, function(x){runGo(x, "BP")})

saveRDS(go, paste0("data/processed/pathway_genesets/go_", model.type,"_any_p.RDS"))


