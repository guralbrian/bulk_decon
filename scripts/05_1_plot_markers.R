libs <- c("tidyverse", "Seurat", "SeuratDisk", "tidyseurat", "viridis", "patchwork", "gt") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

markers <- read.csv("data/processed/single_cell/cluster_markers.csv")
sn.annot <- LoadH5Seurat("data/processed/single_cell/celltype_labeled.h5seurat")

dup.genes <- markers[duplicated(markers$gene),"gene"]

# get highest pvalue genes
genes <- markers |> 
  filter(!(gene %in% dup.genes)) |> 
  group_by(celltype) |>
  arrange(p.value) |> 
  slice_head(n = 6) |> 
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
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
    text = element_text(size = 11)
  ) +
  labs(color = "Average\nExpression", size = "Proportion of\nnuclei expressing")

# Save 
if(!dir.exists("results/5_findMarkers")){
  dir.create("results/5_findMarkers")
}

png(file = "results/5_findMarkers/marker_specificity.png",
    width = 8, 
    height = 3.7,
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
  ggtitle(paste0(top.genes[[i]], " - ", unique(markers$celltype)[[i]])) + NoLegend() + NoAxes()
}

# Save plot to results 
png(file = "results/5_findMarkers/featurePlots.png",
    width = 9.8, 
    height = 6.3,
    units = "in",
    res = 600)

wrap_plots(p.feat)

dev.off()

# Save supp table 
if(!dir.exists("results/supp_figs/")){
  dir.create("results/supp_figs/")
}
options(digits = 3)
tab <- markers |> 
gt(rowname_col = "gene",
   groupname_col = "celltype") |>
  cols_label(
    p.value = "p-value",
    FDR = "False discovery rate",
    summary.logFC = "log-fold-change"
  ) |> 
  cols_align(
    align = "center",
    columns = c(p.value, FDR, summary.logFC)
  ) |> 
  tab_header(
    title = md("Cardiac cell type markers"),
    subtitle = "Top 15 per cell type"
  ) |>  
  opt_row_striping() |> 
  cols_width(gene ~ 90) |>  
  tab_style(
    style = list(
      align = "center",
      cell_fill("grey"),
      cell_text(color = "black", weight = "bold")),
    locations = cells_row_groups())

file = "results/supp_figs/cell_markers"
gtsave(tab, paste0(file, ".html"))

webshot::webshot(url = paste0(file, ".html"), file = paste0(file, ".png"), 
                 vwidth = 550, vheight = n.terms*110, zoom = 3)

