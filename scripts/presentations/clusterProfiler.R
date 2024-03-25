# ClusterProfiler analysis
libs <- c("tidyverse", "clusterProfiler", "DESeq2", 
          "org.Mm.eg.db", "viridis", "stringr", "patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Apply the function
#go <- readRDS("data/processed/pathway_genesets/goadjusted_005.RDS")
go <- readRDS("data/processed/pathway_genesets/gounadjusted_005.RDS")
#go <- go[c(1,2,5,3,4)]
go.joined <- lapply(go, clusterProfiler::simplify, cutoff = 0.6)

plotGO <- function(x, title){
df <- x |> 
  as.data.frame() |> 
  mutate(qscore = -log(p.adjust, base=10),
         desc.wrap = str_to_title(Description) |> 
                    str_wrap(width = 18) |> 
                     factor()) |> 
  arrange(desc(qscore)) |> 
  slice_head(n = 7) 
df$desc.wrap <- factor(df$desc.wrap, levels = rev(df$desc.wrap))

p.ego <- ggplot(df, aes(x = desc.wrap, y = qscore)) +
    geom_bar(stat = "identity", color = "black", fill = "#ed9209") +
    coord_flip() +
    scale_fill_viridis(direction = -1) +
    theme(
      axis.text.y = element_text(hjust = 0.5, size = 12),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "darkgray"),
      legend.position = "none",
      #legend.key.size = unit(1, "in"),
      title = element_text(size = 18)
    ) +
    labs(fill = "Adjusted\np-value",
         y = "Q Score",
         title = title)
}

# Make title lists
titles <- c()
titles[["adjusted"]] <- c("MI", "cmAKO", "MI x cmAKO", "Fibroblasts", "Cardiomyocytes")
titles[["unadjusted"]] <- c("MI", "cmAKO", "MI x cmAKO")

#p.go <- lapply(1:length(go.joined), function(n){plotGO(go.joined[[n]], titles[["adjusted"]][[n]])})
p.go <- lapply(1:length(go.joined), function(n){plotGO(go.joined[[n]], titles[["unadjusted"]][[n]])})

# Make the directory to populate the results in
if(!dir.exists("results/mohlke")){
  dir.create("results/mohlke")
}

# Save plot to results 
png(file = paste0("results/mohlke/unadj_joined_interaction_clusters.png"),
    width = 14, 
    height = 8,
    units = "in",
    res = 600)

wrap_plots(p.go)

dev.off()

