# Plots for Christoph
# 06/22/2023

sn.clust <- LoadH5Seurat( "jensen/data/processed/single_cell/06282023_NoDoublets")

genes <- row.names(sn.clust) |>
  as.data.frame()

plot.dim <- do_DimPlot(sn.clust, reduction = "umap", pt.size = 2)  + 
  NoLegend() + 
  scale_color_brewer(palette = "Dark2")

plots <- c()


plots[["Pln"]] <- sn.clust |>
  do_FeaturePlot("Pln", reduction = "umap", slot = "scale.data", order = T,
                 legend.position = "bottom", legend.width = 3, legend.length = 40) +
  labs(title = "Pln") +
  theme(title = element_text(size = 35),
        legend.title = element_blank() )


plots[["Myh7"]] <- sn.clust |>
  do_FeaturePlot("Myh7", reduction = "umap", slot = "data", order = T,
                 legend.position = "bottom", legend.width = 3, legend.length = 40) +
  labs(title = "Myh7") +
  theme(title = element_text(size = 35),
        legend.title = element_blank() )


plots[["Myh6"]] <- sn.clust |>
  do_FeaturePlot("Myh6", reduction = "umap", slot = "scale.data", order = T,
                 legend.position = "bottom", legend.width = 3, legend.length = 40) +
  labs(title = "Myh6") +
  theme(title = element_text(size = 35),
        legend.title = element_blank() )


plots[["Col1a1"]] <- sn.clust |>
    do_FeaturePlot("Col1a1", reduction = "umap", slot = "scale.data", order = T,
                  legend.position = "bottom", legend.width = 3, legend.length = 40) +
    labs(title = "Col1a1") +
    theme(title = element_text(size = 35),
          legend.title = element_blank() ) 

wrap_plots(plots, 
           ncol = 2)

sn.clust |>
  do_FeaturePlot("Serpina3n", reduction = "umap", slot = "scale.data", order = T,
                 legend.position = "right", legend.width = 3, legend.length = 40) +
  labs(title = "Serpina3n") +
  theme(title = element_text(size = 35),
        legend.title = element_blank())

sn.clust |>
  do_FeaturePlot("H1f0", reduction = "umap", slot = "scale.data", order = T,
                 legend.position = "right", legend.width = 4) +
  labs(title = "H1f0") +
  theme(title = element_text(size = 35),
        legend.title = element_blank())


plot.dim <- do_DimPlot(sn.clust, reduction = "umap", pt.size = 2)  

color_scheme <- brewer.pal(brewer.pal.info["Dark2", "maxcolors"], "Dark2")


library(ggplot2)

# Define the new names for the groups
new_names <- c("Endothelial Cells", "Macrophages", "Fibroblasts", "B Cells", "Cardiomyocytes")

# Create a named vector where names are the old names and the values are the new names
names_vector <- setNames(new_names, levels(sn.clust))

p <- VlnPlot(sn.clust, features = "Serpina3n", slot = "scale.data", 
             idents = levels(sn.clust), cols = color_scheme, sort = "increasing") +
  scale_x_discrete(labels = names_vector) + # This line is for renaming x-axis labels
  labs(
    title = "Serpina3n Expression by Cell Type Cluster",
    x = "Cell Types",
    y = "Scaled Expression"
  ) +
  theme(axis.text.x = element_text(color = "black", size = 25,angle = 0, hjust = 0.5), # Increase size here for larger axis text
      axis.text.y = element_text(color = "black", size = 15), # Add this line for y-axis text size
      axis.ticks.x = element_line(color = "black", size = 1), # Move legend to upper right
      panel.background = element_blank(),
      title = element_text(size = 30),
      axis.title.x =  element_text(vjust= 0.2, size = 30),
      axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
      axis.line.y = element_line(color="black", size = 0.5),
      plot.margin = margin(b = 20)) +
  NoLegend()
# Print the plot

p



# Serpina3n expression corr

matrix <- sn.clust@assays$RNA@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Serpina3n",])
correlations <- apply(matrix_mod,1,function(x){cor(gene,x)})
top.corr.df <- data.frame(genes = names(top.corr),
                          correlations = correlations)

top.corr <- correlations[order(correlations, decreasing = T)] |>
  as.data.frame()
colnames(top.corr) <- "correlation"

# just fibroblasts
sn.fib <- subset(sn.clust, idents = "fibroblasts")
matrix <- sn.fib@assays$RNA@counts
matrix_mod<- as.matrix(t(matrix))
gene<-as.numeric(matrix_mod["Serpina3n",])
correlations.fib <- apply(matrix_mod,1,function(x){cor(gene,x)})
top.corr.fib <- correlations.fib[order(correlations.fib, decreasing = T)] |>
  as.data.frame()

colnames(top.corr.fib) <- "correlation"

write.csv(top.corr, "jensen/results/coexpression/serpina3n/correlations_all")
write.csv(top.corr.fib, "jensen/results/coexpression/serpina3n/correlations_fib")



# potential network method
library("scLink")

count = readRDS(system.file("extdata", "example.rds", package = "scLink"))
genes = readRDS(system.file("extdata", "genes.rds", package = "scLink"))
count.norm = sclink_norm(count, scale.factor = 1e6, filter.genes = FALSE, gene.names = genes)


sn.fib <- subset(sn.clust, idents = "fibroblasts")
matrix <- sn.fib@assays$RNA@counts
matrix_mod<- as.matrix(t(matrix))

count.norm = sclink_norm(matrix_mod, scale.factor = 1e6, filter.genes = T, n = 500 )

networks = sclink_net(expr = count.norm, ncores = 2, lda = seq(0.5, 0.1, -0.05))

