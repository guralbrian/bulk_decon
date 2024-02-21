libs <- c("tidyverse", "Seurat", "SeuratDisk", "tidyseurat") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

markers <- read.csv("data/processed/single_cell/all_markers.csv")
sn <- LoadH5Seurat("data/processed/single_cell/celltype_labeled.h5seurat")


genes <- markers |> 
  group_by(annotation) |> 
  arrange(desc(summary.logFC)) |> 
  slice_head(n = 3) |> 
  pull(gene) 

# Use colourblind-friendly colours
friendly_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC")

# Set theme
my_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
      )
  )

test <- sn |> 
  join_features(features = genes, shape = "long") |> 
  select(c(.cell, .feature, .abundance_RNA, cell.type)) |> 
  mutate(expressed = case_when(
    .abundance_RNA == 0 ~ "no",
    .default = "yes"
  )) |> 
  group_by(cell.type, .feature) 
# Need to get % of cell expressing the gene and average expression of those cells