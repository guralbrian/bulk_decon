# load libraries 
library(dplyr) # load first since we'll need dplyr::lapply()
libs <- c("Seurat", "DESeq2", "patchwork","SeuratDisk", "MuSiC", "reshape2", 
          "tidyverse", "SingleCellExperiment","harmony") # list libraries here
lapply(libs, require, character.only = T)

source("scripts/functions/decon_all.R") # load our own functions

##### Load data #####
gtex.sn <- LoadH5Seurat("data/processed/internal/sn_gtex_no_de.h5seurat")
gtex.bk <- read.csv("data/processed/internal/gtex_lv_counts_summed.csv", check.names = F, row.names = 2)
gtex.bk <- gtex.bk[,-1]


# cluster and remove mito-heavy nuclei
sn.clust <- ClusterSeurat(gtex.sn, 
                              res = 0.2,
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
dim.ft  <- c("Broad.cell.type", "seurat_clusters" )
feat.ft <- c("PercentMito", "PercentRibo")

plotUMAP(data = sn.clust,
         dim.ft = dim.ft,
         feat.ft = feat.ft)

# select bad clusters/cells
bad.clust <- mitoProps(sn.clust,
                  cutoff = 0.01)
print(bad.clust)

bad.cells <- WhichCells(sn.clust, 
                        idents = names(bad.clust[bad.clust > 20])) # sets a 20% cutoff

sn.clust <- subset(gtex.sn, cells = bad.cells, invert = T)


# redo QC with removed cells 
sn.clust <- ClusterSeurat(sn.clust, 
                          res = 0.2,
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

plotUMAP(data = sn.clust,
         dim.ft = dim.ft,
         feat.ft = feat.ft)

#### Music #####

# music_prop needs a SCE format & ExpressionSet  
#Idents(sn.clust) <- sn.clust$Broad.cell.type
sn.sce <- as.SingleCellExperiment(sn.clust, assay = "RNA")

# music_prop needs exprs matrix
bulk.es <- ExpressionSet(assayData = as.matrix(gtex.bk))
bulk.es <- exprs(bulk.es)
# estimate bulk composition with music
decon <- music_prop(bulk.mtx = bulk.es, sc.sce = sn.sce, markers = NULL,
                    clusters = "seurat_clusters", samples = "Participant.ID")

# turn music output into graph-friendly dataframe
# also adds ratios of nuclei clusters from subject-matched samples
props <- processProps(decon.obj = decon,
                           seurat.obj = sn.clust)
decon.melt = melt(props, id.vars = c("names", "group"))
colnames(decon.melt) = c('Sub',"Type", 'CellType', 'Prop')
decon.melt$CellType = factor(decon.melt$CellType, levels = unique(decon.melt$CellType))
decon.melt$Prop <- as.numeric(decon.melt$Prop)
##### Visualize #####

matches <- levels(sn.clust$Participant.ID)
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
    scale_x_discrete(labels = matches)

  
aspect_ratio <- 1/0.4375

#ggsave(filename = "results/music_gtex_matched", height = 7 , width = 7 * aspect_ratio,
#       device = "png")


# all gtex, ordered by cm content

rank <- props[order(props[which(props$group == "Bulk"),4]),length(props)-1]


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
  xlab("GTEx Individuals")
