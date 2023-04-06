# load libraries 
library(dplyr) # load first since we'll need dplyr::lapply()
libs <- c("Seurat", "DESeq2", "patchwork","SeuratDisk", "MuSiC", "reshape2", "tidyverse", "SingleCellExperiment","harmony") # list libraries here
lapply(libs, require, character.only = T)

source("scripts/functions/decon_all.R") # load our own functions

##### Load data #####
gtex.sn <- LoadH5Seurat("data/processed/internal/sn_gtex_lv_match.h5seurat")
good.genes <- read.csv("data/processed/internal/gtex_stable_genes.csv")
gtex.bk <- read.csv("data/processed/internal/gtex_lv_counts_summed.csv", check.names = F, row.names = 1)

# Make new seurat with stable genes 
counts <- GetAssayData(gtex.sn, slot="counts", assay="RNA")[good.genes$x,] |>
            CreateSeuratObject()

# add meta data
meta.features <- colnames(gtex.sn@meta.data)
for(i in meta.features){
  counts <- AddMetaData(counts, gtex.sn@meta.data[[i]], col.name = i)
}

# reassign and clear junk
gtex.sn <- counts
rm(counts)
gc()

# cluster and remove mito-heavy nuclei
sn_gtex_mito <- ClusterSeurat(gtex.sn, 
                              res = 0.05,
                              subset = T,
                              min.rna.ft = 200,
                              max.rna.ft = 2500,
                              min.rna.ct = 800,
                              max.mt.pt = 0.05,
                              max.rb.pt = 0.05,
                              scrublet_score  = 0.4,
                              nfeatures = 2000,
                              harmony = T,
                              regress.by = "batch")

# standard visuals for clusters and mitocondrial contaminants
dim.ft  <- c("Participant.ID","seurat_clusters", "Broad.cell.type")
feat.ft <- c("PercentMito", "PercentRibo", "exon_ratio")

dim.plots  <- vector("list", length(dim.ft))
feat.plots <- vector("list", length(feat.ft))

for(i in seq_along(dim.plots)){
 dim.plots[[i]]<- DimPlot(sn_gtex_mito, reduction = 'umap', group.by = dim.ft[[i]])
}


for(i in seq_along(feat.plots)){
 feat.plots[[i]] <-  FeaturePlot(sn_gtex_mito, 
                        reduction = "umap", 
                        features = feat.ft[[i]],
                        pt.size = 0.4, 
                        order = TRUE,
                        label = TRUE)  
}

patchwork::wrap_plots(c(dim.plots, feat.plots))


#### select bad cluster and redo QC without nuclei within it ####
median.mito <- sn_gtex_mito@meta.data |>
  group_by(seurat_clusters) |>
  summarize(medMito = median(PercentMito)) |>
  as.data.frame()

mean.mito <- sn_gtex_mito@meta.data |>
  group_by(seurat_clusters) |>
  summarize(medMito = mean(PercentMito)) |>
  as.data.frame()


bad.clusters <- as.numeric(mean.mito[as.numeric(mean.mito$medMito) > 0.005, 1])-1
bad.cells <- WhichCells(sn_gtex_mito, idents = bad.clusters)

sn_gtex_trim <- subset(sn_gtex_new, cells = bad.cells, invert = T)


# redo QC with removed cells 
sn_gtex_trim <- ClusterSeurat(sn_gtex_trim, 
                              res = 0.05,
                              subset = T,
                              percentMito = 0.005,
                              harmony = F,
                              featuresRNAmin = 200,
                              countRNA = 500)

DimPlot(sn_gtex_trim, reduction = 'umap', group.by = "seurat_clusters")
FeaturePlot(sn_gtex_trim, 
            reduction = "umap", 
            features = "PercentMito",
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

VlnPlot(sn_gtex_trim, features = "PercentMito", split.by = "seurat_clusters")

#### Music #####

# music_prop needs a SCE format  
sn.sce <- as.SingleCellExperiment(sn_gtex_trim, assay = "RNA")

# This is a test to see how simple log2 norm of counts works

#counts(sn.sce) <- log2(counts(sn.sce)+1)
#bk_gtex <- log2(bk_gtex + 1)

# music_prop needs exprs matrix
bulk.es <-ExpressionSet(assayData = as.matrix(bk_gtex))
bulk.es.trim <- exprs(bulk.es)[,-1] 

#only include bulk rna seq which matches sn
sn_p <- c("13N11","15RIE", "1ICG6")
matches <- unique(grep(paste(sn_p,collapse="|"), 
                       colnames(bulk.es.trim), value=TRUE))
bulk.es.trim <- bulk.es.trim[,matches] 

# estimate bulk composition with music
clusters.keep <- c(0:9)
decon <- music_prop(bulk.mtx = bulk.es.trim, sc.sce = sn.sce, markers = NULL,
                    clusters = "seurat_clusters", samples = "Participant.ID",
                    select.ct = clusters.keep)
# turn music output into graph-friendly dataframe
# also adds ratios of nuclei clusters from subject-matched samples
props <- processProps(decon.obj = decon,
                           seurat.obj = sn_gtex_mito,
                           subjects = sn_pid)
decon.melt = melt(props)
colnames(decon.melt) = c('Sub',"Type", 'CellType', 'Prop')
decon.melt$CellType = factor(decon.melt$CellType, levels = unique(decon.melt$CellType))

##### Visualize #####

matches <- unique(grep(paste(sn_p,collapse="|"), 
                       decon.melt$Sub, value=TRUE))
legend_labels <- new.cluster.ids[(1+as.numeric(colnames(props)[1:6]))]

# plot comparing subject-matched bulk and sn estimates

ggplot(decon.melt[which(decon.melt$Sub %in% matches),], 
       aes(x=Sub, y=Prop, fill=CellType))  +
    geom_bar(stat='identity', 
             position = "fill", 
             width = 0.5,
             color = "black")+
    scale_fill_discrete(name = "Cell Type")+
    ylab("Proportion") +
    facet_grid(~Type, scales = "free_x") +
    theme(
      axis.text.x = element_text(color = "black", size = 15),
      strip.text = element_text(size = 20),
      title = element_text(size = 20),
      legend.text = element_text(size = 18),
      axis.ticks.x = element_blank()) +
    xlab("Samples")+
    scale_x_discrete(labels = c(sn_pid))

  
aspect_ratio <- 1/0.4375

#ggsave(filename = "results/music_gtex_matched", height = 7 , width = 7 * aspect_ratio,
#       device = "png")


# all gtex, ordered by cm content

rank <- props[order(props[,1]),7]


ggplot(decon.melt[which(decon.melt$Type == "Bulk"),], 
       aes(x=factor(Sub, levels = rank), y=Prop, fill=CellType))  +
  geom_bar(stat='identity', 
           position = "fill", 
           width = 1)+
  scale_fill_discrete(name = "Cell Type")+
  ylab("Proportion") +
  facet_grid(~Type, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    strip.text = element_blank(),
    title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.ticks.x = element_blank()) +
  xlab("GTEx Individuals")+
  scale_x_discrete(labels = c(sn_pid))

aspect_ratio <- 1/0.4375

#ggsave(filename = "results/music_gtex_all", height = 7 , width = 7 * aspect_ratio,
#      device = "png")


# Cluster Markers code

cluster2.markers <- FindAllMarkers(sn_gtex_trim,
                                   logfc.threshold = 0.583,
                                   min.pct = 0.25, 
                                   min.diff.pct = 0.25,
                                   only.pos = TRUE)


clusters <- cluster2.markers |>
  group_by(cluster) |>
  slice_max(n = 200, order_by = avg_log2FC)
genes(0)

new.cluster.ids <- c("Fibroblast II",
                     "Endothelium",
                     "Macrophage",
                     "Ventricular CM I",
                     "Pericyte",
                     "Adipocyte")

sn_gtex_new <- RenameIdents(sn_gtex_trim, new.cluster.ids)



##### looking at gene weights in deconvolution ####

genes.8 <- genes(8)
gene.weight <- as.data.frame(decon$Weight.gene)
gene.weight$mean <- rowMeans(gene.weight)
gene.weight.sub <- gene.weight[which(rownames(gene.weight) %in% genes.8$gene),]
decon$Weight.gene[rownames(decon$Weight.gene) %in% genes(8),]


FeaturePlot(sn_gtex_mito, 
            reduction = "umap", 
            features = "FREM1",
            pt.size = 0.4, 
            order = TRUE,
            label = TRUE)
