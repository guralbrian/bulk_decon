libs <- c("tidyverse", "Seurat", "SeuratDisk", "tidyseurat", "viridis") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

markers <- read.csv("data/processed/single_cell/all_markers.csv")
sn.annot <- LoadH5Seurat("data/processed/single_cell/celltype_labeled.h5seurat")

# get highest pvalue genes
genes <- markers |> 
  group_by(annotation) |> 
  arrange(p.value) |> 
  slice_head(n = 3) |> 
  arrange(annotation) |> 
  pull(gene) 

# get df of cells with gene info
sn.markers <- sn.annot |> 
  join_features(features = genes, shape = "long") |> 
  select(c(.cell, .feature, .abundance_RNA, cell.type)) |> 
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
  geom_point(aes(size = cells_expressing + cells_expressing*0.5)) +
  geom_point(aes(size = cells_expressing, color = cells_expressing)) +
  scale_color_viridis(option="plasma", direction = -1) +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2)
  ) +
  labs(color = "Average\nExpression", size = "Proportion of\nnuclei expressing")

# Save 
png(file = "results/5_findMarkers/marker_specificity.png",
    width = 10, 
    height = 5,
    units = "in",
    res = 500)

p.mark

dev.off()

