# Load libraries

libs <- c("Seurat", "ggplot2", "DESeq2", "patchwork","SeuratDisk", "MuSiC", "reshape2",
          "tidyverse", "SingleCellExperiment", "harmony", "SCpubr", 
          "AUCell", "viridis", "gplots", "scales", "ggrepel", "gridExtra", "scCustomize",
          "httr","matrixStats", "scran", "scuttle", "scater", "DropletUtils", "scDblFinder") # list libraries here

lapply(libs, require, character.only = T)
source("jensen/scripts/functions/decon_all.R")

# load individual seurats
doublet.path <- "jensen/data/processed/single_cell/no_doublets/"
datasets <- list.files(doublet.path)

sn.list <- lapply(datasets, function(x){LoadH5Seurat(paste0(doublet.path, x))})

# filter each of them by quantile values
sn.list <- lapply(sn.list, function(x){FilterByQuantile(x)})
# Join into one

sn.all <- merge(sn.list[[1]], sn.list[[2]]) |>
            merge(sn.list[[3]]) |>
            merge(sn.list[[4]])
sn.clust <- sn.all |>
    subset(isDoublet == "no") |>
    ClusterSeurat(subset = F, res = 0.03, regress.by = c("origin"))

# standard visuals for clusters and mitocondrial contaminants
dim.ft  <- c("seurat_clusters")
feat.ft <- "PercentMito"
plotUMAP(data = sn.clust,
         dim.ft = dim.ft,
         feat.ft = feat.ft,
         nrow = 1, 
         ncol = 2)

do_DimPlot(sn.clust, group.by = "isDoublet", reduction = "umap", label = T, repel = T) + NoLegend()
do_FeaturePlot(sn.clust, features = "DoubletScore", reduction = "umap") + NoLegend()


VlnPlot(sn.clust, "DoubletScore", group.by = "orig.ident")


# Annotate to cell types

markers.broad <- read.csv("jensen/data/processed/external/mclellan_2020/mclellan_cell_markers_broad.csv")

# Split 'markers' into separate data frames for each unique subcluster
subclusters <- split(markers.broad, markers.broad$cluster)
# Extract the gene names from each subcluster data frame
subcluster_genes <- lapply(subclusters, function(x) x$gene)
# Create a list with named elements corresponding to each subcluster
geneSets <- setNames(subcluster_genes, names(subclusters))
# make expression matrix of single cell
my.expr <-  sn.clust@assays$RNA@counts
# Find  AUC for each cell by each type
cells_AUC <- AUCell_run(my.expr, geneSets)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE)

# get the thresholds and report when/how they're passdd
selectedThresholds <- getThresholdSelected(cells_assignment)

auc.val <- cells_AUC@assays@data@listData$AUC |>
  as.data.frame() |>
  t()

# Apply process_cell_id function to each row of auc.val using lapply()
broad_cell_types_list <- lapply(1:nrow(auc.val), function(i) processCellId(auc.val[i, ]))

# Convert the list to a data frame
broad.cell.types <- data.frame(BroadCellType = unlist(broad_cell_types_list), stringsAsFactors = FALSE)
rownames(broad.cell.types) <- rownames(auc.val)

sn.clust <- AddMetaData(sn.clust, broad.cell.types)



#### Visualize Annotation corr ####
# Create a contingency table
contingency_table <- table(c(sn.clust$BroadCellType),
                           c(sn.clust$seurat_clusters))

# Perform the chi-squared test
chi_squared_test <- chisq.test(contingency_table)

# Calculate the standardized residuals
observed_frequencies <- contingency_table
expected_frequencies <- chi_squared_test$expected
standardized_residuals <- (observed_frequencies - expected_frequencies) / sqrt(expected_frequencies)


col_palette <- viridis(25)
# Create a heatmap of the standardized residuals
heatmap.2(standardized_residuals, 
          trace = "none", 
          col = col_palette,
          main = "Residuals, All Cell Types",
          xlab = "Seurat Clusters",
          ylab = "All Cell Types",
          margins = c(8, 12),
          key.title = "Standardized Residuals",
          srtCol = 0,
          cexRow = 1.2)

# assign labels which are max and above threshold
threshold   <- 80

sn.clust <- AssignAndFilterClusters(sn.clust)


# bulk
bulk <- read.csv("jensen/data/processed/jensen_rau_froese_cpm", row.names = 1)
bulk_pheno <- read.csv("jensen/data/processed/jensen_rau_froese_pheno", row.names = 1)
bulk.es <- ExpressionSet(assayData = as.matrix(bulk[,colnames(bulk) != "wt_Lx3"]))
bulk.es <- exprs(bulk.es)
sn.list <- c()
sn.list[["original_sn"]] <- LoadH5Seurat("jensen/data/processed/single_cell/06152023.h5seurat")
sn.list[["new_sn"]] <- AssignAndFilterClusters(sn.clust)
props <- lapply(sn.list, function(x){EstimateCellTypeProportions(x, bulk.es, for.aitchison = T)}) 

#### Aitchison ####


library("coda.base")

# rows = samples , cols = proportions
sim.samples <- bulk_pheno[which(bulk_pheno$fraction != "whole_tissue") ,c(1,3,4)]

# assume that each sample differs in purity
# cms = 0.95, ecs_rau = 0.7, ecs_froese = 0.9, fbs_rau = 0.7, fbs_froese = 0.7


cell.dict <- data.frame(bulk.pheno = unique(sim.samples$fraction),
                        cell.types = c("cardiomyocytes", "fibroblasts", "endothelial cells"))

sim.samples$expected.purity <- 0.95

sim.samples[which(sim.samples$fraction == "Fibroblasts"),4] <- 0.9
sim.samples[which(sim.samples$fraction == "Endothelial Cells"),4] <- 0.9
sim.samples[which(sim.samples$fraction == "Endothelial Cells" & sim.samples$origin == "Rau"),4] <- 0.7

# make dataframes with expected compositions for each fraction (70 - 95% purity)

sim.fractions <- lapply(sn.list, function(sn){CreateSimFractions(sn, sim.samples, cell.dict, purity.adjustment = 1)})


# Generate 20 values sampled from a normal distribution between 0.9 and 1.05
set.seed(123) # For reproducibility
purity_adjustments <- rnorm(20, mean = 0.95, sd = 0.05)
purity_adjustments <- ifelse(purity_adjustments < 0.9, 0.9, purity_adjustments)
purity_adjustments <- ifelse(purity_adjustments > 1.05, 1.05, purity_adjustments)

# Use lapply to iterate over sn.combos and purity_adjustments
sim.fractions <- lapply(sn.list, function(sn) {
  lapply(purity_adjustments, function(pa) {
    CreateSimFractions(sn, sim.samples, cell.dict, purity.adjustment = pa)
  })
})

aitch_vals <- lapply(seq_along(sim.fractions), function(i) {
  lapply(seq_along(sim.fractions[[i]]), function(j) {
    CalculateAitchisonDistance(sim.fractions[[i]][[j]], props[[i]])
  })
})

names(aitch_vals) <- c("original_sn", "new_sn")
cell.names <- bulk_pheno[which(bulk_pheno$fraction != "whole_tissue"),]
#### Plotting ####

# Create a new dataframe from aitch_vals
df <- do.call(rbind, lapply(seq_along(aitch_vals), function(i){
  do.call(rbind, lapply(seq_along(aitch_vals[[i]]), function(j){
    if(nrow(aitch_vals[[i]][[j]]) > 0){
      temp_df <- cbind(data.frame(type=names(aitch_vals)[i], 
                                  aitch_vals[[i]][[j]], check.names = FALSE),
                       id = rownames(aitch_vals[[i]][[j]]))
      temp_df
    } else{
      NULL
    }
  }))
}))

# Merge with cell.names
df <- merge(df, cell.names, by = "id", all.x = TRUE)

# Convert to a factor for ordered plotting
df$type <- reorder(df$type.x, -df$aitchison, FUN = sum)
df$aitchison.inverse <- df$aitchison

# Plot
ggplot(df, aes(x = type, y = aitchison)) +
  geom_boxplot(aes(fill = fraction), 
               width = 0.5, color = "black", outlier.shape = NA) +
  geom_point(aes(shape = fraction, color = origin),
             position = position_jitterdodge(jitter.width = 0.05, jitter.height = 0, dodge.width = 0.5), size = 2, alpha = 0.5) +
  theme(axis.text.x = element_text(color = "black", size = 15, angle = 30, vjust = 0.6),
        axis.text.y = element_text(color = "black", size = 10),
        axis.ticks.x = element_line(color = "black", size = 1),
        panel.background = element_blank(),
        strip.text = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(color="black", size = 0.5)) +
  labs(x = "Cell Type", y = "Aitchison Value", fill = "Cell Fraction", color = "Origin", shape = "Cell Fraction") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(16, 17, 18)) # use different shapes for each fraction


plotUMAP(sn.clust,
         feat.ft = c("Serpina3n", "Col1a1"))
sn.clust |>
  subset(origin == "tabula_muris") |>
do_FeaturePlot( "Serpina3n", reduction = "umap") + 
  labs(title = "Serpina3n expression in Tabula Muris")

do_DimPlot(sn.clust, reduction = "umap", label = T, repel = T) + NoLegend()


VlnPlot(sn.clust, c("Serpina3n"), group.by = "origin")

row.names(sn.clust)
