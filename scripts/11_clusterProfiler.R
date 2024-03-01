# ClusterProfiler analysis
libs <- c("tidyverse", "clusterProfiler", "DESeq2", "org.Mm.eg.db") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)


# Load the results and expression matrix
# Results from model with cell type proportions
#  "treatment_TAC_vs_Sham" "genotype_KO_vs_WT""clr.Cardiomyocytes" "clr.Fibroblast" "treatmentTAC.genotypeKO"
adj.res <- readRDS("data/processed/models/adjusted_de_interaction.RDS")  |> 
  results(name="genotype_KO_vs_WT") |> 
  as.data.frame()

# Results from model without cell type proportions
unadj.res <- readRDS("data/processed/models/unadjusted_de_interaction.RDS")  |> 
  results(name="treatmentTAC.genotypeKO") |> 
  as.data.frame()

# Function to perform GO enrichment with clusterProfiler
runGo <- function(data, onto){
# Get a list of significant genes
sig.genes <- data |> 
  filter(padj < 0.05) |> 
  row.names()

# Convert common gene names to ENSEMBLE IDs for clusterProfiler
gene.df <- bitr(sig.genes, fromType = "SYMBOL",
                toType = c("ENSEMBL"),
                OrgDb = org.Mm.eg.db)

# Check for enrichment within biological process gene clusters
ego <- enrichGO(gene          = gene.df$ENSEMBL,
                OrgDb         = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont           = onto,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
}

# Apply the function
go.adj <- runGo(adj.res, "BP")
#go.unadj <- runGo(unadj.res, "BP")


mutate(go.adj, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
