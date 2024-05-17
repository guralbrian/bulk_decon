libs <- c("gridExtra", "tidyverse", "Seurat", "SeuratDisk", "tidyseurat", "viridis", "RColorBrewer","patchwork", "cowplot") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

sn.annot <- LoadH5Seurat("data/processed/single_cell/celltype_labeled.h5seurat")


# get df of cells with gene info
sn.markers <- sn.annot |> 
  join_features(features = "Cryab", shape = "long") |> 
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
  geom_point(aes(size = cells_expressing *1.3)) +
  geom_point(aes(size = cells_expressing, color = cells_expressing)) +
  scale_color_viridis(option = "magma") +
  theme(
    axis.text.y = element_text(angle = 30, size = 10),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2)
  ) +
  labs(color = "Average\nExpression", size = "% nuclei expressing")

p.mark
# Make feature plots
# loop through markers and plot for each
p.feat <- FeaturePlot(sn.annot, "Cryab", order = T, pt.size = 1.2) + 
    ggtitle("Cryab expression") + NoLegend()

p.feat

# Color scheme
cols <-  brewer.pal(length(levels(sn.annot$cell.type)), "Set2")
names(cols) <- levels(sn.annot$cell.type)

p.umap <- DimPlot(sn.annot, group.by = "cell.type", cols = cols, pt.size = 1) |>
  LabelClusters(id = "cell.type", size = 6, repel = T, box = T, force = 100) +
  ggtitle(NULL) +
  theme_cowplot(font_size = 18) +
  NoLegend()

p.vil <- VlnPlot(sn.annot, group.by = "cell.type", features = "Cryab", 
                 sort = T, cols = cols) +
                NoLegend() +
                ggtitle(NULL)  +
                labs(x = NULL,
                     y = "Cryab expression") +
              theme(
                axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5, size = 15),
                axis.title.y = element_text(size = 16)
              ) 
p.vil


# Save plot to results 
png(file = "results/misc/tijia_sn.png",
    width = 6, 
    height = 12,
    units = "in",
    res = 600)


grid.arrange(p.umap, p.feat, p.vil,
             heights = c(2, 2, 2),
             layout_matrix = rbind(c(1),
                                   c(2),
                                   c(3)))


dev.off()

