# ClusterProfiler analysis
libs <- c("tidyverse", "clusterProfiler", "DESeq2", 
          "org.Mm.eg.db", "viridis", "stringr", "patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Get commandArgs
args <- commandArgs(trailingOnly = TRUE)
model.type <-  as.character(args[1])

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
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE) #|> 
        #clusterProfiler::simplify(cutoff = 0.6)
}

# Apply the function
go <- lapply(res, function(x){runGo(x, "BP")})

saveRDS(go, paste0("data/processed/pathway_genesets/go", model.type,"_005.RDS"))

plotGO <- function(x, title){
df <- x |> 
  as.data.frame() |> 
  mutate(qscore = -log(p.adjust, base=10),
         desc.wrap = str_to_title(Description) |> 
                    str_wrap(width = 18) |> 
                     factor()) |> 
  arrange(desc(qscore)) |> 
  slice_head(n = 10) 
df$desc.wrap <- factor(df$desc.wrap, levels = rev(df$desc.wrap))

p.ego <- ggplot(df, aes(x = desc.wrap, y = qscore, fill = p.adjust)) +
    geom_bar(stat = "identity", color = "black") +
    coord_flip() +
    scale_fill_viridis(direction = -1) +
    theme(
      axis.text.y = element_text(hjust = 0.5, size = 12),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "darkgray"),
      legend.position = "bottom"
    ) +
    labs(fill = "Adjusted\np-value",
         y = "Q Score",
         title = title)
}

# Make title lists
titles <- c()
titles[["adjusted"]] <- c("MI", "cmAKO", "Fibroblasts", "Cardiomyocytes", "MI x cmAKO")
titles[["unadjusted"]] <- c("MI", "cmAKO", "MI x cmAKO")

p.go <- lapply(1:length(go), function(n){plotGO(go[[n]], titles[[model.type]][[n]])})

# Make the directory to populate the results in
if(!dir.exists("results/11_clusterProfiler")){
  dir.create("results/11_clusterProfiler")
}

# Save plot to results 
png(file = paste0("results/11_clusterProfiler/", model.type, "_interaction_clusters.png"),
    width = 21, 
    height = 18,
    units = "in",
    res = 600)

wrap_plots(p.go)

dev.off()

