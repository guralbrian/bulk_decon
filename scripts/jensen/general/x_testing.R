# scCATCH #####
library("scCATCH")
obj <- createscCATCH(data = sn.clust@assays$RNA@counts, 
                     cluster = as.character(sn.clust$seurat_clusters))

tissues  <- cellmatch[which(cellmatch$species == "Mouse"),2] |>
              unique() |>
              as.data.frame()
# make marker
cellmatch_new <- cellmatch[cellmatch$species == "Mouse" & cellmatch$tissue %in% c("Heart", "Heart muscle"), ]
obj <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = cellmatch_new)
obj <- findcelltype(obj)

levels(obj@celltype$cell_type)

#### EasyCellType ####
library(org.Mm.eg.db)
library(AnnotationDbi)

markers <- FindAllMarkers(sn.clust, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)

markers$entrezid <- mapIds(org.Mm.eg.db,
                           keys=markers$gene, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
markers <- na.omit(markers)

markers_sort <- data.frame(gene=markers$entrezid, cluster=markers$cluster, 
                           score=markers$avg_log2FC) %>% 
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 
input.d <- as.data.frame(markers_sort[, 1:3])

annot.GSEA <- easyct(input.d, db="cellmarker", species="Mouse", 
                     tissue= c("Heart", "Heart muscle"), p_cut=0.8,
                     test="GSEA", scoretype = "std")

plot_dot(test="GSEA", data = annot.GSEA)
cell.type <- data.frame(type1 = rep(NA, length(levels(sn.clust))),
                        type2 = rep(NA, length(levels(sn.clust))),
                        type3 = rep(NA, length(levels(sn.clust))),
                        p1    = rep(NA, length(levels(sn.clust))),
                        p2    = rep(NA, length(levels(sn.clust))),
                        p3    = rep(NA, length(levels(sn.clust))))
count <- 0
annots <- rep(NA, length(levels(sn.clust)))
for(i in levels(sn.clust)){
  count <- count + 1
  layers <- length(annot.GSEA[[as.character(i)]]$ID)
  layers.use <- min(c(layers, 3))
  cell.type[count, c(1:layers.use)] <- annot.GSEA[[as.character(i)]]$ID[c(1:layers.use)]
  cell.type[count, c(4:(layers.use+3))] <- annot.GSEA[[as.character(i)]]$pvalue[c(1:layers.use)]
  if(layers == 1){
    annots[[count]] <- cell.type$type1[[count]]
  }else{
  if(cell.type$p1[[count]]<cell.type$p2[[count]]*0.5 ){
    annots[[count]] <- cell.type$type1[[count]]
  }else{
    annots[[count]] <- paste(cell.type[count,c(1:layers)], collapse =  "/")}}
}

annots <- make.unique(annots)







#### SCINA  ####

library(SCINA)


#### AUCell ####

#! you left off here. You were reading in a list of markers for mouse heart cells found in McLellan 2022. 
# load in the markers and use them in AUCell to annotate your seurat clusters
# the name of the file needs to be changed, you ref the wrong paper (doi.org/10.1161/CIRCULATIONAHA.119.045115 is right)

library(AUCell)
# load markers
markers <- read.csv("data/processed/external/skelly_2020/mclellen_cell_markers.csv")
markers.broad <- read.csv("data/processed/external/skelly_2020/mclellen_cell_markers_broad.csv")
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
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
# make a cell assignment table
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
# get the thresholds and report when/how they're passdd
selectedThresholds <- getThresholdSelected(cells_assignment)

cellAUCdf <- matrix(NA, 
                    nrow = length(colnames(cells_AUC)), 
                    ncol = length(rownames(cells_AUC))) |>
                as.data.frame()
rownames(cellAUCdf) <- colnames(cells_AUC)
colnames(cellAUCdf) <- rownames(cells_AUC)
# get boolean values for each cell type assignment and add the metadata
for(i in row.names(cells_AUC)) {
  print(i)
cellAUCdf[,i] <- t(getAUC(cells_AUC[i,]) >  selectedThresholds[i])
sn.clust <- AddMetaData(sn.clust, cellAUCdf[,i] , col.name = i)
}

# now I want to write a part that looks at cells with >1 tag and picks a single cell type
# this will be by looking at the threshold for each type and seeing how it's own AUC compares to it as a ratio
# i.e. if cell_1 has an AUC of 0.3 for CMs and 0.4 for ECs but the thresholds 
#      are 0.1 and 0.3, it'll be tagged as a CMs since 0.3/0.1 > 0.4/0.3

auc.val <- cells_AUC@assays@data@listData$AUC |>
              as.data.frame() |>
              t()


# Create an empty data frame to store the results
max_ratios_cell_types <- data.frame(CellType = character(), stringsAsFactors = FALSE)

# Iterate through each row (cell ID) of the auc.val matrix
for (i in 1:nrow(auc.val)) {
  # Get the cell ID row
  row_values <- auc.val[i, ]
  # Filter the columns with values greater than 0
  filtered_values <- row_values[row_values > 0]
  # Calculate the ratio between the value in auc.val and the corresponding cell type's threshold value from selectedThresholds
  ratios_for_cell_id_list <- lapply(names(filtered_values), function(cell_type) {
    filtered_values[cell_type] / selectedThresholds[[cell_type]]
  })
  ratios_for_cell_id <- unlist(ratios_for_cell_id_list, use.names = TRUE)
  names(ratios_for_cell_id) <- names(filtered_values)
  # Filter the ratios above 1
  filtered_ratios <- ratios_for_cell_id[ratios_for_cell_id > 1]
  # Find the cell type corresponding to the largest ratio for the current cell ID
  if (length(filtered_ratios) > 0) {
    max_ratio_cell_type <- names(filtered_ratios)[which.max(filtered_ratios)]
  } else {
    max_ratio_cell_type <- NA
  }
  # Store the corresponding cell type
  max_ratios_cell_types <- rbind(max_ratios_cell_types, data.frame(CellType = max_ratio_cell_type, stringsAsFactors = FALSE))
}

# Set the rownames to the cell IDs and remove the CellID column
rownames(max_ratios_cell_types) <- rownames(auc.val)

sn.clust <- AddMetaData(sn.clust, max_ratios_cell_types)


# cramers V

library(vcd)
# Assuming your data is stored in a data frame called sn.clust with columns BroadCellType and seurat_clusters
# Create a contingency table
# Assuming your data is stored in a data frame called sn.clust with columns BroadCellType and seurat_clusters
# Create a contingency table
contingency_table <- table(sn.clust$BroadCellType, sn.clust$seurat_clusters)

# Create a mosaic plot
mosaic(contingency_table, shade = TRUE, legend = TRUE, split_vertical = T,main = "Mosaic Plot of BroadCellType and Seurat Clusters")


# Assuming your data is stored in a data frame called sn.clust with columns BroadCellType and seurat_clusters
# Create a contingency table
contingency_table <- table(sn.clust$BroadCellType, sn.clust$seurat_clusters)

library(ggcorrplot)
model.matrix(~0+., data=data.frame(sn.clust$BroadCellType, sn.clust$seurat_clusters)) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2)
