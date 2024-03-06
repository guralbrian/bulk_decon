# ClusterProfiler analysis
libs <- c("tidyverse", "clusterProfiler", "DESeq2", 
          "org.Mm.eg.db", "viridis", "stringr", "patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)


# Load the results and expression matrices
adj.res <- readRDS("data/processed/models/adjusted_de_interaction.RDS")  
adj.names <- resultsNames(adj.res)[-1]

# make a list of each result
adj.res <- lapply(adj.names, function(x){
  results(adj.res, name=x) |> 
    as.data.frame()
})

unadj.res <- readRDS("data/processed/models/unadjusted_de_interaction.RDS")  
unadj.names <- resultsNames(unadj.res)[-1]

unadj.res <- lapply(unadj.names, function(x){
  results(unadj.res, name=x) |> 
    as.data.frame()
})


# Function to perform GO enrichment with clusterProfiler
runGo <- function(data, onto){
# Get a list of significant genes
sig.genes <- data |> 
  filter(padj < 0.05 ) |> #& abs(log2FoldChange) >= 0.263) |> 
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
                readable      = TRUE) |> 
        clusterProfiler::simplify(cutoff = 0.6)
}

# Apply the function
go.adj   <- lapply(adj.res,   function(x){runGo(x, "BP")})
go.unadj <- lapply(unadj.res, function(x){runGo(x, "BP")})

plotGO <- function(x, title){
df <- x |> 
  as.data.frame() |> 
  mutate(qscore = -log(p.adjust, base=10),
         desc.wrap = str_to_title(Description) |> 
                    str_wrap(width = 18) |> 
                     factor()) |> 
  arrange(desc(qscore)) |> 
  slice_head(n = 6) 
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
adj.titles <- c("CAD", "cmAKO", "Cardiomyocytes", "Fibroblasts", "CAD x cmAKO")
unadj.titles <- c("CAD", "cmAKO", "CAD x cmAKO")

p.adj <- lapply(1:length(go.adj), function(n){plotGO(go.adj[[n]], adj.titles[[n]])})
p.unadj <- lapply(1:length(go.unadj), function(n){plotGO(go.unadj[[n]], unadj.titles[[n]])})

# Save plot to results 
png(file = "results/11_clusterProfiler/adj_interaction_clusters.png",
    width = 21, 
    height = 14,
    units = "in",
    res = 300)

wrap_plots(p.adj)

dev.off()

png(file = "results/11_clusterProfiler/unadj_interaction_clusters.png",
    width = 21, 
    height = 7,
    units = "in",
    res = 300)

wrap_plots(p.unadj)

dev.off()

