libs <- c("tidyverse", "Seurat", "SeuratDisk", "tidyseurat", "viridis", "patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

markers <- read.csv("data/processed/single_cell/cluster_markers.csv")
sn.annot <- LoadH5Seurat("data/processed/single_cell/celltype_labeled.h5seurat")

dup.genes <- markers[duplicated(markers$gene),"gene"]
# get highest pvalue genes
genes <- markers |> 
  filter(!(gene %in% dup.genes)) |> 
  group_by(celltype) |> 
  arrange(desc(summary.logFC)) |> 
  slice_head(n = 3) |> 
  arrange(celltype) |> 
  pull(gene) 

# get df of cells with gene info
sn.markers <- sn.annot |> 
  join_features(features = genes, shape = "long") |> 
  dplyr::select(c(.cell, .feature, .abundance_RNA, cell.type)) |> 
  mutate(expressed = case_when(
    .abundance_RNA == 0 ~ 0,
    .default = 1
  )) |> 
  group_by(cell.type, .feature) |> 
  summarize(cells_expressing = sum(expressed)/n(), 
            ave_exp = mean(.abundance_RNA))

# Set factor level orders
sn.markers$.feature <- factor(sn.markers$.feature, levels=genes)

cell.list <- levels(Idents(sn.annot))[order(levels(Idents(sn.annot)))] |> rev()
sn.markers$cell.type <- factor(sn.markers$cell.type, levels=cell.list)

# get % of cell expressing the gene and average expression of clusters
p.mark <- sn.markers |> 
  ggplot(aes(x = .feature, y = cell.type)) +
  geom_point(aes(size = cells_expressing * 1.2)) +
  geom_point(aes(size = cells_expressing, color = cells_expressing)) +
  scale_color_viridis(option = "magma") +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2)
  ) +
  labs(color = "Average\nExpression", size = "Proportion of\nnuclei expressing")

# Save 
if(!dir.exists("results/5_findMarkers")){
  dir.create("results/5_findMarkers")
}

png(file = "results/5_findMarkers/marker_specificity.png",
    width = 8, 
    height = 4.5,
    units = "in",
    res = 500)

p.mark

dev.off()

# Make feature plots
top.genes <- markers |> 
  group_by(celltype) |> 
  arrange(desc(summary.logFC)) |> 
  slice_head(n=1) |> 
  pull(gene)

# loop through markers and plot for each
p.feat <- c()
for(i in 1:length(top.genes)) {
p.feat[[i]] <- FeaturePlot(sn.annot, top.genes[[i]]) + 
  ggtitle(paste0(top.genes[[i]], " - ", unique(markers$celltype)[[i]])) + NoLegend()
}

# Save plot to results 
png(file = "results/5_findMarkers/featurePlots.png",
    width = 9.8, 
    height = 6.3,
    units = "in",
    res = 600)

wrap_plots(p.feat)

dev.off()

