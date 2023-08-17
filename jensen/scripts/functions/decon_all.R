#' Filter a Seurat object by stably 
#' expressed genes between bulk and 
#' single nucleus counts data.
#' 
#' @param seurat.obj       Seurat object containing the single nucleus data.
#' @param bulk.obj         Counts matrix for the bulk data. Format by samples as
#'                         columns and genes as rows. Matches genes by row
#'                         names/genes.
#' @param min.rna.features Minimum number of RNA features to filter by.
#' @param max.rna.features Maximum number of RNA features to filter by.
#' @param min.rna.count    Minimum read count to filter RNA by.
#' @param max.mt.percent   Maximum mitochondrial gene percentage of total single
#'                         nucleus reads to filter cells by.
#' @param max.dims         The number of PCA dimensions to use for generating 
#'                         UMAP reduction.
#' @param visualize        If TRUE, make DimPlot() for the final, filtered 
#'                         single nucleus data. 
#' @param lfc.threshold    Log2 fold-change to consider a gene differentially
#'                         expressed when comparing pseudobulk and real bulk.
#' @param resolution       Sets resolution for clustering nuclei after removing
#'                         DE genes.                       
#' @returns A numeric vector.
#' 
#' @export
FilterBulkSingleNucleus <- function(seurat.obj, 
                                    bulk.obj,
                                    min.rna.ft = 200,
                                    max.rna.ft = 2500,
                                    min.rna.ct = 800,
                                    max.mt.pt  = 0.05,
                                    group      = "orig.ident",
                                    lfc.thresh = 0.583,
                                    change     = "greaterAbs") {
  # seurat qc 
  seurat.obj <- subset(seurat.obj, 
                       subset = nFeature_RNA   > min.rna.ft     & 
                                nFeature_RNA   < max.rna.ft     &
                                nCount_RNA     > min.rna.ct     & 
                                scrublet_score < scrublet.score &
                                PercentRibo   <= max.rb.pt      &
                                PercentMito   <= max.mt.pt)
  
  pseudo.bulk        <- AggregateExpression(seurat.obj, 
                                            group.by = group,
                                            assays = 'RNA',
                                            slot = "counts", 
                                            return.seurat = FALSE)
  pseudo.bulk      <- as.data.frame(pseudo.bulk$RNA)
  colnames(pseudo.bulk) <- paste0(colnames(pseudo.bulk), "-pb")
  
  # match sn and bulk transcripts
  bulk.sn <- merge(pseudo.bulk, bulk.obj, by.x = "row.names", by.y = "row.names")
  genes <- bulk.sn[,1]
  bulk.sn <-  sapply(bulk.sn[,c(2:length(bulk.sn))], as.numeric)
  rownames(bulk.sn) <- genes
  
  #make meta data
  meta.data <- data.frame("id" = colnames(bulk.sn),
                          "type" = as.factor(c(
                            rep("single.nucleus",length(pseudo.bulk)),
                            rep("bulk",length(bulk.obj)))))
  # create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = bulk.sn,
                                colData = meta.data, 
                                design = ~type) 
  
  # diff expression analysis
  dds <- DESeq(dds)
  result_de <- results(dds, lfcThreshold = lfc.thresh, altHypothesis = change,
                       contrast = c("type", "bulk", "single.nucleus"))
  good.genes <- subset(result_de, padj > 0.05)
  good.genes <- rownames(good.genes)

  return(good.genes)
}

ClusterSeurat <- function(seurat.obj,
                          subset = T,
                          min.rna.ft = 200,
                          max.rna.ft = 2500,
                          min.rna.ct = 800,
                          max.mt.pt  = 0.05,
                          regress.by = "orig.ident",
                          harmony    = T,
                          res        = 0.2,
                          nfeatures  = 2000,
                          drop.levels = FALSE){
  set.seed(100)
  if(subset == T){
    seurat.obj <- subset(seurat.obj, 
                         subset = nFeature_RNA   > min.rna.ft     & 
                                  nFeature_RNA   < max.rna.ft     &
                                  nCount_RNA     > min.rna.ct     &
                                  PercentMito   <= max.mt.pt)
  }
  seurat.obj <- seurat.obj |>
    NormalizeData(verbose = F) |>
    FindVariableFeatures(verbose = F, nfeatures = nfeatures) |>
    ScaleData(verbose = F) |>
    RunPCA(verbose = F) 
  
  if(harmony == T){
    if(drop.levels ==T){
    seurat.obj@meta.data[[regress.by]] <- droplevels(seurat.obj@meta.data[[regress.by]]) 
    }
    seurat.obj <- RunHarmony(seurat.obj, 
                             group.by.vars = regress.by,
                             verbose = F,
                             project.dim = F)
  }
  
  # find elbow
    # Determine percent of variation associated with each PC
    pct <- seurat.obj[["pca"]]@stdev / sum(seurat.obj[["pca"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    pcs <- min(co1, co2)
  
  if(harmony == T){
    seurat.obj <- seurat.obj |>
      FindNeighbors(dims = 1:pcs, reduction = "harmony", verbose = F) |>
      FindClusters(resolution = res, verbose = F) |>
      RunUMAP(dims = 1:pcs, reduction = "harmony", verbose = F)
    return(seurat.obj)
  }else{
    seurat.obj <- seurat.obj |>
      FindNeighbors(dims = 1:pcs, reduction = "pca", verbose = F) |>
      FindClusters(resolution = res, verbose = F) |>
      RunUMAP(dims = 1:pcs, reduction = "pca", verbose = F)
    return(seurat.obj)
  }}

# a function to return the top lfc genes from 
genes <- function(x){
  print(as.data.frame(clusters[which(clusters$cluster == x), 7]), row.names = F, print.keys = F)
}

processProps <- function(decon.obj,
                         seurat.obj,
                         subjects,
                         subject.slot = "Participant.ID"){
  decon.melt = melt(decon.obj$Est.prop.weighted)
  colnames(decon.melt) = c('Sub', 'CellType', 'Prop')
  decon.melt$CellType = factor(decon.melt$CellType, levels = unique(decon.melt$CellType))
  
  ##### Add SN ground truths back in ####
  props <- as.data.frame(decon.obj[[1]])
  props <- props[,as.character(seq(0,length(props)-1))] # reorders columns
  props$names <- rownames(props) 
  active.idents <- as.data.frame(table(seurat.obj@active.ident, seurat.obj@meta.data[[subject.slot]]))
  
  # exclude sn cell types not used in decon
  pops.decon <- colnames(decon.obj$Var.prop)
  pops.sn <- levels(active.idents$Var1)
  active.idents <- active.idents[active.idents$Var1 %in% pops.decon,]
  active.idents <- active.idents  |>
    group_by(Var2)|>
    mutate(proportion = Freq/sum(Freq))
  for(i in unique(active.idents$Var2)){
    props <- rbind(props, c(as.vector(t(active.idents[active.idents$Var2 == i, 4])), i))
  }
  
  rownames(props) <- c(rownames(decon.obj$Est.prop.weighted),unique(active.idents$Var2))
  props$group <- c(rep("Bulk",length(rownames(decon.obj$Est.prop.weighted))),
                   rep("Single Nucleus",length(unique(active.idents$Var2))))
  
  return(props)
}

optimMusic <- function(param,
                       i = 1,
                       seurat.obj,
                       bulk.obj,
                       return.props = F,
                       return.seurat = T,
                       subset.samples = F){
  # expand the test parameters 
  test.params <- param[i,]
  params.names <- colnames(param)
  
  for(j in 1:length(params.names)){
    assign(params.names[[j]], test.params[1,j])
  }
  
  # Do DE analysis if applicable 
  if(remove.genes == T){
    if(subset.samples == T){
      matches <- unique(grep(paste(sn_p,collapse="|"), 
                             colnames(bulk.obj), value=TRUE))
      good.genes <- FilterBulkSingleNucleus(seurat.obj = seurat.obj, 
                                            bulk.obj = bulk.obj[,matches], 
                                            max.mt.percent = mt.percent1,
                                            group = "Participant.ID",
                                            lfc.threshold = log2(fold.change))
    }else{
      good.genes <- FilterBulkSingleNucleus(seurat.obj = seurat.obj, 
                                            bulk.obj = bulk.obj, 
                                            max.mt.percent = mt.percent1,
                                            group = "Participant.ID",
                                            lfc.threshold = log2(fold.change))
    }
  }
  # remake seurat by excluding DE genes
  counts <- GetAssayData(seurat.obj, slot="counts", assay="RNA")   
  
  if(remove.genes == T){
    genes.filter <- counts[(rownames(counts) %in% good.genes), ]
  }else{
    genes.filter <- counts
  }
  
  sn_gtex_new <- CreateSeuratObject(counts = genes.filter)
  
  meta.features <- c("Participant.ID",
                     "nFeature_RNA",
                     "nCount_RNA",
                     "PercentMito")
  
  for(i in meta.features){
    sn_gtex_new <- AddMetaData(sn_gtex_new,seurat.obj@meta.data[[i]], col.name = i)
  }
  
  sn_gtex_new@meta.data$Participant.ID <- droplevels(sn_gtex_new@meta.data$Participant.ID)
  
  ##### cluster and remove mito-heavy nuclei #####
  sn_gtex_mito <- ClusterSeurat(sn_gtex_new, 
                                res = res1,
                                subset = T,
                                percentMito = mt.percent2,
                                harmony = harmony,
                                featuresRNAmin = 200,
                                countRNA = 500)
  if(remove.clusters == T) {
    # remove clusters and redo clustering
    mean.mito <- sn_gtex_mito@meta.data |>
      group_by(seurat_clusters) |>
      summarize(medMito = mean(PercentMito)) |>
      as.data.frame()
    
    good.clusters <- as.numeric(mean.mito[as.numeric(mean.mito$medMito) < 0.003, 1])-1
    good.cells <- WhichCells(sn_gtex_mito, idents = good.clusters)
    sn_gtex_trim <- subset(sn_gtex_new, cells = good.cells, invert = F)
    
    
    # redo QC with removed cells 
    sn_gtex_trim <- ClusterSeurat(sn_gtex_trim, 
                                  res = res2,
                                  subset = T,
                                  percentMito = mt.percent3,
                                  harmony = harmony,
                                  featuresRNAmin = 200,
                                  countRNA = 500)
    
    
    sn.sce <- as.SingleCellExperiment(sn_gtex_trim, assay = "RNA")
  }else{
    sn.sce <- as.SingleCellExperiment(sn_gtex_mito, assay = "RNA")
    sn_gtex_trim <- sn_gtex_mito
  }
  if(return.seurat == T){
    return(sn_gtex_trim)
  }
  # music_prop needs exprs matrix
  bulk.es <-ExpressionSet(assayData = as.matrix(bulk.obj)) |>
    exprs()
  
  # estimate bulk composition with music
  decon <- music_prop(bulk.mtx = bulk.es, sc.sce = sn.sce, markers = NULL,
                      clusters = "seurat_clusters", samples = "Participant.ID")
  
  # turn music output into graph-friendly dataframe
  # also adds ratios of nuclei clusters from subject-matched samples
  props <- processProps(decon.obj = decon,
                        seurat.obj = sn_gtex_trim,
                        subjects = sn_pid)
  if(return.props == T){
    return(props)
  }
  decon.melt <- melt(props, id.vars = c("names", "group"))
  colnames(decon.melt) = c('Sub',"Type", 'CellType', 'Prop')
  decon.melt$CellType = factor(decon.melt$CellType, levels = unique(decon.melt$CellType))
  
  decon.1 <- decon.melt[which(decon.melt$Type == "Bulk"),c(1,3,4)]
  colnames(decon.1) <- c("names", "variable", "bulk.est")
  
  decon.2 <- decon.melt[which(decon.melt$Type == "Single Nucleus"),c(1,3,4)]
  colnames(decon.2) <- c("names", "variable", "sn.est")
  decon.2$names <- sapply(strsplit(decon.2$names, "-"), "[[", 2)
  decon.all <- merge(decon.1,decon.2,by=c("names","variable"))
  decon.all$bulk.est <- as.numeric(decon.all$bulk.est)
  decon.all$sn.est <- as.numeric(decon.all$sn.est)
  # make a function to  get ratios of cell types
  good.rows <- apply(decon.all[,c(3,4)], 1, function(row) all(row !=0 ))
  decon.sub <- decon.all[good.rows,]
  ratios.bulk <- data.frame()
  clusters <- as.numeric(levels(unique(decon.sub$variable)))
  cell.combos <- expand.grid(type1 = clusters,
                             type2 = clusters)
  cell.combos <- cell.combos[which(cell.combos$type1 != cell.combos$type2),]
  for(x in unique(decon.all$names)){
    data <- decon.all[which(decon.all$names == x),]
    ratios <- data.frame(ratio = NA,
                         type1 = NA,
                         type2 = NA)
    for(y in 1:nrow(cell.combos)){
      ratios[y,1] <- as.numeric(data[which(data$variable == cell.combos[y,1]),3]) / 
        as.numeric(data[which(data$variable == cell.combos[y,2]),3])
      ratios[y,c(2,3)] <- cell.combos[y,c(1,2)]
    }
    ratios$sub <- rep(x, nrow(ratios))
    ratios.bulk <- rbind(ratios.bulk, ratios)
  }
  
  # ratios of sn cell proportions 
  ratios.sn <- data.frame()
  for(x in unique(decon.all$names)){
    data <- decon.all[which(decon.all$names == x),]
    ratios <- data.frame(ratio = NA,
                         type1 = NA,
                         type2 = NA)
    for(y in 1:nrow(cell.combos)){
      ratios[y,1] <- as.numeric(data[which(data$variable == cell.combos[y,1]),4]) / 
        as.numeric(data[which(data$variable == cell.combos[y,2]),4])
      ratios[y,c(2,3)] <- cell.combos[y,c(1,2)]
    }
    ratios$sub <- rep(x, nrow(ratios))
    ratios.sn <- rbind(ratios.sn, ratios)
  }
  
  
  ratios.bulk$ratio.sn <- ratios.sn$ratio
  # get rid of infinite and NA values
  bad.rows <- is.infinite(ratios.bulk$ratio) | is.nan(ratios.bulk$ratio) | 
    is.infinite(ratios.bulk$ratio.sn) | is.nan(ratios.bulk$ratio.sn)
  ratios.bulk <- ratios.bulk[!(bad.rows),]
  
  # correlations
  r.ratios.spearman <- cor(ratios.bulk$ratio,ratios.bulk$ratio.sn, method = "spearman")^2
  r.ratios.pearson <- cor(ratios.bulk$ratio,ratios.bulk$ratio.sn, method = "pearson")^2
  rmse <- decon.all |>
    summarise(rmse = sqrt(mean((as.numeric(bulk.est) - as.numeric(sn.est))^2)))
  r.single.spearman <- cor(as.numeric(decon.all[,3]),as.numeric(decon.all[,4]), method = "spearman")^2
  r.single.pearson <- cor(as.numeric(decon.all[,3]),as.numeric(decon.all[,4]), method = "pearson")^2
  results <- c(rmse[1,1],
               r.ratios.spearman,
               r.ratios.pearson,
               r.single.spearman,
               r.single.pearson,
               test.params) |>
    t()
  return(results)
}

plotUMAP <- function(data,
                     dim.ft = NULL,
                     feat.ft = NULL,
                     feat.labels =NULL,
                     vln.ft = NULL,
                     vln.groups = NULL,
                     width = NULL,
                     height = NULL,
                     ncol = NULL,
                     nrow = NULL,
                     design = NULL,
                     label.size = 4,
                     font.size = NULL,
                     name = NULL) {
  dim.plots  <- vector("list", length(dim.ft))
  feat.plots <- vector("list", length(feat.ft))
  vln.plots <- vector("list", length(vln.ft)*length(vln.groups))
  print(name)
  for(i in seq_along(dim.plots)){
    dim.plots[[i]]<- do_DimPlot(data, 
                             reduction = 'umap', 
                             group.by = dim.ft[[i]],
                             label = T,
                             repel = T,
                             font.size = font.size,
                             label.size = label.size) + 
                          NoLegend() + 
                          ggtitle(name) +
                          theme(plot.title = element_text(size = font.size*2, hjust = 0.5))
  } 
  
  for(i in seq_along(feat.plots)){
    feat.plots[[i]] <-  SCpubr::do_FeaturePlot(data, 
                                    reduction = "umap", 
                                    features = feat.ft[[i]],
                                    pt.size = 0.4, 
                                    order = TRUE,
                                    font.size = font.size,
                                    label.size = label.size) + 
                                    ggtitle(if(
                                      !is.null(feat.labels[[i]])){
                                        paste0("*",feat.ft[[i]],"*", ", ", feat.labels[[i]])
                                           }else{
                                        feat.ft[[i]]}) +
                                    theme(plot.title = ggtext::element_markdown())
  }
  for(i in seq_along(vln.ft)){
    for(j in seq_along(vln.groups)){
    vln.plots[[i*j]] <-  SCpubr::do_ViolinPlot(data, 
                                               features = vln.ft[[i]], 
                                               group.by = vln.groups[[j]],
                                               font.size = font.size) + 
                                               ggtitle(paste0(vln.ft[[i]], " by ", vln.groups[[j]]))
    }
  }
  return(patchwork::wrap_plots(c(dim.plots, feat.plots, vln.plots), 
                               widths = width,
                               heights = height,
                               ncol = ncol,
                               nrow = nrow,
                               design = design))
}

# returns the proportion of nuclei in each cluster above a value of meta data 
mitoProps <- function(data,
                      cutoff = 0.03,
                      origins.considered = "rau"){
                      sn.mito <- subset(data, subset = PercentMito > cutoff &
                                          origin %in% origins.considered)
                      sn.mito.tb <- table(sn.mito$seurat_clusters)
                      sn.origin <- subset(data, subset = origin %in% origins.considered)
                      sn.tb <- table(sn.origin$seurat_clusters)
                      clust.pcs <- (sn.mito.tb / sn.tb) * 100 |>
                        round(digits = 1)
                      return(clust.pcs)
}

# Create a function to process a single row (cell ID) from auc.val
processCellId <- function(row_values) {
  row_values <- auc.val[1,]
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
  return(max_ratio_cell_type)
}

ClusterByMeta <- function(seurat.obj,
                          sub.group,
                          subset = T,
                          min.rna.ft = NULL,
                          max.rna.ft = NULL,
                          min.rna.ct = NULL,
                          max.mt.pt  = NULL,
                          meta.groups = "orig.ident",
                          regress.by = "orig.ident",
                          harmony    = T,
                          res        = 0.2,
                          nfeatures  = 2000,
                          drop.levels = FALSE,
                          doublet_detection = T,
                          ambient_correction = T){
  set.seed(100)
  if(subset == T){
    # Get unique origins
    seurat.sub <- seurat.obj[
      which(
        seurat.obj@meta.data[[meta.group]] == sub.group),]
    
    # Calculate default filtering values if not provided
    if(is.null(min.rna.ft)) min.rna.ft <- quantile(seurat.sub$nFeature_RNA, 0.05)
    if(is.null(max.rna.ft)) max.rna.ft <- quantile(seurat.sub$nFeature_RNA, 0.95)
    if(is.null(min.rna.ct)) min.rna.ct <- quantile(seurat.sub$nCount_RNA, 0.05)
    if(is.null(max.mt.pt)) max.mt.pt <- quantile(seurat.sub$PercentMito, 0.95)
    
    # Subset the Seurat object based on the calculated thresholds
    seurat.sub <- subset(seurat.sub, 
                         subset = nFeature_RNA   > min.rna.ft     & 
                           nFeature_RNA   < max.rna.ft     &
                           nCount_RNA     > min.rna.ct     &
                           PercentMito   <= max.mt.pt)
    # Normalize
    set.seed(100)
    seurat.sub <- seurat.sub |>
      NormalizeData(verbose = F) |>
      FindVariableFeatures(verbose = F, nfeatures = nfeatures) |>
      ScaleData(verbose = F) |>
      RunPCA(verbose = F) 
    print("Normalized")
    # find elbow
    # Determine percent of variation associated with each PC
    pct <- seurat.sub[["pca"]]@stdev / sum(seurat.sub[["pca"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    pcs <- min(co1, co2)
    
    # cluster
    set.seed(100)
    seurat.sub <- seurat.sub |>
      FindNeighbors(dims = 1:pcs, reduction = "pca", verbose = F) |>
      FindClusters(resolution = res, verbose = F) |>
      RunUMAP(dims = 1:pcs, reduction = "pca", verbose = F)
    print("U = MAP'd")
    # Doublet Finder Preprocessing
    
    # pK identification (no ground-truth)
    sweep.list <- paramSweep_v3(seurat.sub, PCs = 1:pcs)
    sweep.stats <- summarizeSweep(sweep.list)
    bcmvn <- find.pK(sweep.stats)
    print("pK = found")
    
    # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
    bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
    optimal.pk <- bcmvn.max$pK
    optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
    
    ## Homotypic doublet proportion estimate
    annotations <- seurat.sub@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations) 
    nExp.poi <- round(optimal.pk * nrow(seurat.sub@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
    # run DoubletFinder
    seurat.sub <- doubletFinder_v3(seu = seurat.sub, 
                                   PCs = 1:pcs, 
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
  }
}


RemoveDoublets <- function(sub_list, seurat, directory){
  
  sn.sub <- subset(x = seurat, subset =  orig.ident == sub_list)
  
  
  # Find empty droplets 
  # make sce
  sce.sub <- as.SingleCellExperiment(sn.sub)
  
  bcrank <- barcodeRanks(counts(sce.sub))
  uniq <- !duplicated(bcrank$rank)
  
  # Now at each plot, you wrap the plot call with png(), and dev.off() to save the image
  plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
  abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
  abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
  legend("bottomleft", legend=c("Inflection", "Knee"), col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
  
  # Save plot
  ggsave(filename = paste0(directory,"/droplets_elbow_", sub_list, ".png"), units = "px",dpi=300, width = 1000, height = 900 )
  print(paste0("saved elbowplot of ", sub_list))
  
  # set cutoff limit and consider droplets below it to be empty
  set.seed(100)
  limit <- metadata(bcrank)$inflection   
  e.out <- emptyDrops(counts(sce.sub), lower=limit, test.ambient=TRUE)
  
  sce2 <- sce.sub[,which(e.out$FDR <= 0.001)]
  
  
  # Define and remove ambient RNA 
  clusters <- quickCluster(sce2)
  sce2 <- computeSumFactors(sce2, cluster=clusters)
  
  # evaluate ambinat RNA contamination in the empty droplets
  amb <- metadata(e.out)$ambient[,1]
  head(amb)
  
  sce2 <- logNormCounts(sce2)
  set.seed(1000)
  # modeling variables
  dec.pbmc <- modelGeneVarByPoisson(sce2)
  # calcualte top features
  top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)
  #
  set.seed(1000)
  # Evaluate PCs
  sce2 <- denoisePCA(sce2, subset.row=top.pbmc, technical=dec.pbmc)
  # make UMAP plot
  sce2 <- runUMAP(sce2, dimred="PCA")
  g <- buildSNNGraph(sce2, k=25, use.dimred = 'PCA')
  clust <- igraph::cluster_walktrap(g)$membership
  colLabels(sce2) <- factor(clust)
  
  # Save plotUMAP
  
  scater::plotUMAP(sce2, colour_by="label")
  ggsave(filename = paste0(directory, "/UMAP_sce2_", sub_list, ".png"), units = "px", dpi =300, width = 1000, height = 900 )
  
  
  stripped <- sce2[names(amb),]
  out <- removeAmbience(counts(stripped), ambient= amb, groups = colLabels(stripped))
  print(paste0("removed ambience from ", sub_list))
  # load correccted counts into scce object
  counts(stripped, withDimnames=FALSE) <- out
  stripped <- logNormCounts(stripped)
  
  # Find Doublets ####
  
  dbl.dens <- computeDoubletDensity(stripped, #subset.row=top.mam, 
                                    d=ncol(reducedDim(stripped)),subset.row=top.pbmc)
  
  stripped$DoubletScore <- dbl.dens
  
  # Save plotUMAP with DoubletScore
  
  scater::plotUMAP(stripped, colour_by="DoubletScore")
  ggsave(filename = paste0(directory, "/doubletsUMAP_", sub_list, ".png"), units = "px",dpi=300, width = 1000, height = 900 )
  
  # Save plotColData with DoubletScore
  plotColData(stripped, x="label", y="DoubletScore", colour_by="label")+
    geom_hline(yintercept = quantile(colData(stripped)$DoubletScore,0.95),lty="dashed",color="red")
  ggsave(filename = paste0(directory, "/doubletsViolin_", sub_list, ".png"), units = "px",dpi=300, width = 1000, height = 900 )
  
  
  cut_off <- quantile(stripped$DoubletScore,0.95)
  stripped$isDoublet <- c("no","yes")[factor(as.integer(stripped$DoubletScore>=cut_off),levels=c(0,1))]
  
  sn.clean <- as.Seurat(stripped)
  print(paste0("reverted back to seurat: ", sub_list))
  return(sn.clean)
}

FilterByQuantile <- function(seurat.obj,
                             min.rna.ft = NULL,
                             max.rna.ft = NULL,
                             min.rna.ct = NULL,
                             max.mt.pt  = NULL, 
                             pt.remove = 0.1){
  
  seurat.obj$PercentMito <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
  # Calculate default filtering values if not provided
  if(is.null(min.rna.ft)) min.rna.ft <- quantile(seurat.obj$nFeature_RNA, pt.remove)
  if(is.null(max.rna.ft)) max.rna.ft <- quantile(seurat.obj$nFeature_RNA, 1-pt.remove)
  if(is.null(min.rna.ct)) min.rna.ct <- quantile(seurat.obj$nCount_RNA, pt.remove)
  if(is.null(max.mt.pt)) max.mt.pt <- quantile(seurat.obj$PercentMito, 1-pt.remove)
  
  # Subset the Seurat object based on the calculated thresholds
  seurat.obj <- subset(seurat.obj, 
                       subset = nFeature_RNA   > min.rna.ft     & 
                         nFeature_RNA   < max.rna.ft     &
                         nCount_RNA     > min.rna.ct     &
                         PercentMito   <= max.mt.pt)
}


AssignAndFilterClusters <- function(seurat, res.thresh = 0.4, ratio.thresh = 2, min.cell = 400) {
  # Exclude clusters below specific max residual threshold
  max <- apply(standardized_residuals, 2, max)
  max.exclude <- max < res.thresh
  
  # Exclude clusters with few cells 
  bad.clusts <- table(Idents(seurat))[table(Idents(seurat)) < min.cell]
  
  # Get final exclusion list
  all.exclude <- max.exclude 
  max.exclude[is.na(max.exclude)] <- TRUE
  # Find the index of the largest value in each column of the standardized_residuals
  max.indices <- apply(standardized_residuals, 2, which.max)
  
  max.indices[max.exclude] <- NA
  
  # Get the corresponding BroadCellType for each Seurat cluster
  assigned.cell.types <- rownames(standardized_residuals)[max.indices]
  assigned.cell.types[!is.na(assigned.cell.types)] <- make.unique(assigned.cell.types[!is.na(assigned.cell.types)], sep = "_")
  # Rename the clusters
  names(assigned.cell.types) <- names(max.indices)
  seurat <- RenameIdents(seurat, assigned.cell.types)
  
  # Remove very small clusters 
  bad.clusts <- names(table(Idents(seurat)))[table(Idents(seurat)) < min.cell]
  Idents(seurat)[which(Idents(seurat) %in% bad.clusts)] <- NA
  return(seurat)
}

EstimateCellTypeProportions <- function(seurat, bulk.es, for.aitchison = F,sn.individuals, cells_exclude = c("unlabeled", "NA"), bulk.log2 = T, marker = NULL) {
  
  # Convert to SingleCellExperiment
  seurat_sce <- as.SingleCellExperiment(seurat, assay = "RNA")
  
  # Exclude specified clusters
  cells <- levels(Idents(seurat))
  cells <- cells[!(cells %in% cells_exclude) & !is.na(cells)]
  
  # Use MuSiC to estimate cell type proportions
  set.seed(100)
  if(bulk.log2 == T){
    bulk.es <- 2^bulk.es
  }
  decon <- music_prop(bulk.mtx = bulk.es, sc.sce = seurat_sce, markers = marker,
                      clusters = "ident", samples = sn.individuals,
                      select.ct = cells)
  if(for.aitchison == T){
    return(decon)
  }
  # Turn MuSiC output into graph-friendly dataframe
  decon.melt = reshape2::melt(decon$Est.prop.weighted)
  colnames(decon.melt) = c('Sub', 'CellType', 'Prop')
  #decon.melt$combination <- paste(unique(seurat$origin), collapse = "_")
  return(decon.melt)
}

CreateSimFractions <- function(seurat, sim_samples, cell_dict, purity.adjustment = 1, included.cells = NULL) {
  set.seed(100)
  # Get the cell types
  cell_types <- included.cells
  
  # Initialize the sim.fractions dataframe
  sim_fractions <- matrix(nrow = nrow(sim_samples), ncol = length(cell_types)) %>%
    as.data.frame()
  row.names(sim_fractions) <- sim_samples$id
  colnames(sim_fractions) <- cell_types
  
  # Fill in the dataframe
  for(i in 1:nrow(sim_fractions)){
    name1 <- sim_samples[i, 2]
    major <- cell_dict[which(cell_dict$bulk.pheno == name1), 2]
    sim_fractions[i, major] <- sim_samples[i,4] * purity.adjustment
    sim_fractions[i, cell_types[which(cell_types != major)]] <- (1 - (sim_samples[i,4] * purity.adjustment)) / length(cell_types)
  }
  
  return(sim_fractions)
}

# Function to compare estimated and simulated ("real") compositions of cell type fractions
CalculateAitchisonDistance <- function(sim_fractions, est_fractions) {
  ests <- est_fractions$Est.prop.weighted
  ests <- ests[sim.samples$id,]
  ests[ests == 0] <- 0.05 * 0.65
  aitch_vals <- data.frame(aitchison = rep(NA, dim(ests)[1]))
  rownames(aitch_vals) <- rownames(ests)
  cell_types <- colnames(ests)
  
  for(i in rownames(sim_fractions)){
    aitch_vals[i,] <- coda.base::dist(
      rbind(sim_fractions[i,cell_types], ests[i,cell_types]), 
      method = 'aitchison')[1]
  }
  
  return(aitch_vals)
}


MakeContingency <- function(metadata_1 = sn.clust.new$BroadCellType, 
                            metadata_2 = sn.clust.new$seurat_clusters,
                            na.arg = "no"){
  cont_tbl <- table(metadata_1, metadata_2, useNA = na.arg)
  # Perform the chi-squared test
  chi_squared_test <- chisq.test(cont_tbl)
  
  # Calculate the standardized residuals
  obs_freq <- cont_tbl
  exp_freq <- chi_squared_test$expected
  standardized_residuals <- (obs_freq - exp_freq) / sqrt(exp_freq)
  
  return(standardized_residuals)
}


AssignAnnotations <- function(seurat, markers, n_markers = 5, n_cores = 1, n_cells = 1000){
  seurat <- seurat[,sample(names(seurat$seurat_clusters), n_cells)]
  markers.sub <- markers %>%
    group_by(cluster) %>%
    top_n(n_markers, gene) |>
    as.data.frame()
  # Split 'markers' into separate data frames for each unique subcluster
  subclusters <- split(markers.sub, markers.sub$cluster)
  # Extract the gene names from each subcluster data frame
  subcluster_genes <- lapply(subclusters, function(x) x$gene)
  # Create a list with named elements corresponding to each subcluster
  geneSets <- setNames(subcluster_genes, names(subclusters))
  # make expression matrix of single cell
  my.expr <-  seurat@assays$RNA@counts
  # Find  AUC for each cell by each type
  cells_AUC <- AUCell_run(my.expr, geneSets, BPPARAM = if(!is.null(n_cores)){BiocParallel::MulticoreParam(n_cores)})
  
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE, nCores = n_cores)
  
  comments <- lapply(names(cells_assignment), function(cell_type){cells_assignment[[cell_type]]$aucThr$comment})
  good_cells <- names(cells_assignment)[comments == ""]
  # get the thresholds and report when/how they're passdd
  selectedThresholds <- getThresholdSelected(cells_assignment)
  
  auc.val <- cells_AUC@assays@data@listData$AUC |>
    t()|>
    as.data.frame() 
  auc.val <- auc.val[,good_cells]
  auc.val$seurat <- seurat$seurat_clusters[rownames(auc.val)]
  
  return(auc.val)
  
}