# ClusterProfiler analysis
libs <- c("tidyverse", "clusterProfiler", "DESeq2", 
          "org.Mm.eg.db", "viridis", "stringr", "patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load the results and expression matrices
unadj.res <- readRDS("data/processed/pathway_genesets/go_unadj_005_lfc4.RDS")  
adj.res <- readRDS("data/processed/pathway_genesets/go_adj_005_lfc4.RDS")  
palette <- c("#66C2A5","#FC8D62","#8DA0CB")
plotGO <- function(x, descriptors, title){
  df <- x |> 
    as.data.frame() 
  df[which(df$Description =="1-phosphatidylinositol-3-kinase regulator activity"), "Description"] <- "1-phosphatidyl inositol-3-kinase regulator activity"
  df <- df |>   
    subset(Description %in% descriptors) |> 
    mutate(qscore = -log(p.adjust, base=10),
           desc.wrap = str_to_title(Description) |> 
             str_wrap(width = 16, whitespace_only = F) |> 
             factor()) |> 
    arrange(desc(qscore)) 
  df$desc.wrap <- factor(df$desc.wrap, levels = rev(df$desc.wrap))
  
  p.ego <- ggplot(df, aes(x = desc.wrap, y = qscore, fill = desc.wrap)) +
    geom_bar(stat = "identity", color = "black") +
    coord_flip() +
    scale_fill_manual(values = palette) +
    theme(
      axis.text.y = element_text(hjust = 0.5),
      axis.title.y = element_text(size = 30),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "darkgray"),
      legend.position = "none",
      text = element_text(size = 25),
      plot.title = element_blank()
      
    ) +
    labs(y = "Q Score",
         x = title)
}

adj.res[[2]]@result |> pull(Description)

adj.res <- adj.res[c(1,2,5,3,4)]

cad.unadj <- c("oxidative phosphorylation", "extracellular matrix organization", "wound healing")
ko.unadj <- c( "RNA polymerase II-specific DNA-binding transcription factor binding","transcription coregulator binding","histone deacetylase binding" )
int.unadj <- c("interleukin-1 production","muscle cell differentiation" , "regulation of tumor necrosis factor production")

unadj.des <- c(list(cad.unadj),list(ko.unadj),list(int.unadj))


cad.adj <- c("heart contraction", "leukocyte migration", "transmembrane receptor protein serine/threonine kinase signaling pathway")
ko.adj <- c( "circadian regulation of gene expression","protein refolding","regulation of protein catabolic process")
int.unadj <- c("1-phosphatidylinositol-3-kinase regulator activity", "contractile actin filament bundle", "regulation of smooth muscle cell proliferation")
fib.adj <- c("collagen fibril organization", "sarcoplasmic reticulum calcium ion transport", "response to oxidative stress")
cm.adj <- c("cellular respiration","ATP metabolic process", "mitochondrial gene expression" )

adj.des <- c(list(cad.adj),list(ko.adj),list(int.unadj),list(fib.adj),list(cm.adj))
# Make title lists
adj.titles <- c("LAD Ligation", "cmAKO", "Ligation x cmAKO", "Fibroblasts","Cardiomyocytes")
unadj.titles <- c("CAD", "cmAKO", "CAD x cmAKO")


p.unadj <- lapply(1:length(unadj.res), function(n){plotGO(unadj.res[[n]], unadj.des[[n]], unadj.titles[[n]])})
p.adj <- lapply(1:length(adj.res), function(n){plotGO(adj.res[[n]],adj.des[[n]], adj.titles[[n]])})

# Save plot to results 
png(file = "results/mhi/unadj_interaction_clusters.png",
    width = 6, 
    height = 18,
    units = "in",
    res = 400)

wrap_plots(p.unadj, ncol = 1) + 
  plot_annotation('Composition-blind',theme=theme(plot.title=element_text(hjust=0.5, size = 40)))

dev.off()

png(file = "results/mhi/adj_interaction_clusters.png",
    width = 7, 
    height = 30,
    units = "in",
    res = 400)

wrap_plots(p.adj, ncol = 1) + 
  plot_annotation('Cell-types accounted for',theme=theme(plot.title=element_text(hjust=0.5, size = 40)))


dev.off()

