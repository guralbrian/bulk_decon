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
                                    min.rna.features = 200,
                                    max.rna.features = 2500,
                                    min.rna.count = 800,
                                    max.mt.percent = 0.5,
                                    group = "donor",
                                    lfc.threshold = 0.583,
                                    change = "lessAbs") {
  seurat.obj <- subset(seurat.obj, 
                       subset = nFeature_RNA > min.rna.features & nFeature_RNA < max.rna.features &
                         nCount_RNA > min.rna.count & 
                         PercentMito <= max.mt.percent)
  
  pseudo.bulk        <- AggregateExpression(seurat.obj, 
                                            group.by = group,
                                            assays = 'RNA',
                                            slot = "counts", 
                                            return.seurat = FALSE)
  pseudo.bulk      <- as.data.frame(pseudo.bulk$RNA)
  colnames(pseudo.bulk) <- sapply(strsplit(colnames(pseudo.bulk), "-"), "[[", 2)
  # match sn and bulk transcripts
  colnames(bulk.obj) <- sn_p
  bulk.sn <- merge(pseudo.bulk, bulk.obj, by = 0)
  seqs <- bulk.sn$Row.names
  ids <- c(colnames(pseudo.bulk), colnames(bulk.obj))
  # make the counts into integers
  bulk.sn <- bulk.sn[,2:length(bulk.sn)] |>
    lapply(as.integer) |>
    as.data.frame()
  rownames(bulk.sn) <- seqs
  #make meta data
  meta.data <- data.frame("id" = ids,
                          "type" = as.factor(c(
                            rep("single.nucleus", length(pseudo.bulk)),
                            rep("bulk",length(bulk.obj)))))
  # create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = bulk.sn,
                                colData = meta.data, 
                                design = ~type + id) 
  
  # diff expression analysis
  dds <- DESeq(dds)
  result_de <- results(dds, lfcThreshold = lfc.threshold, altHypothesis = "greaterAbs",
                       contrast = c("type", "bulk", "single.nucleus"))
  good.genes <- subset(result_de, padj > 0.05)
  good.genes <- rownames(good.genes)
  # this is the moneymaker step 
  # filter out the transcripts not represented in the res_filter obj
  return(good.genes)
}

ScaleMakePseudobulk <- function(seurat.obj, 
                                min.rna.features = 200,
                                max.rna.features = 2500,
                                min.rna.count = 800,
                                max.mt.percent = 5,
                                group) {
  seurat.obj <- subset(seurat.obj, 
                       subset = nFeature_RNA > min.rna.features & nFeature_RNA < max.rna.features &
                         nCount_RNA > 800 & 
                         percent_mito < 5)
  seurat.obj <- seurat.obj |>
    FindVariableFeatures() |>
    ScaleData() 
  
  pseudo.bulk        <- AggregateExpression(seurat.obj, 
                                            group.by = group,
                                            assays = 'RNA',
                                            slot = "counts", 
                                            return.seurat = FALSE)
  pseudo.bulk      <- as.data.frame(pseudo.bulk$RNA)
}


ClusterSeurat <- function(seurat.obj,
                          res = 0.2,
                          subset = F,
                          percentMito = 0.05,
                          featuresRNAmin = 400,
                          featuresRNAmax = 3000,
                          countRNA = 800,
                          nfeatures = 2000,
                          harmony = F){
  if(subset == T){
    seurat.obj <- subset(seurat.obj, 
                         nFeature_RNA > featuresRNAmin & 
                           nFeature_RNA < featuresRNAmax &
                           nCount_RNA > countRNA &
                           PercentMito <= percentMito)
  }
  seurat.obj <- seurat.obj |>
    NormalizeData(verbose = F) |>
    FindVariableFeatures(verbose = F, nfeatures = nfeatures) |>
    ScaleData(verbose = F) |>
    RunPCA(verbose = F) 
  
  if(harmony == T){
    seurat.obj <- RunHarmony(seurat.obj, 
                             group.by.vars = "Participant.ID")
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
      RunUMAP(dims = 1:pcs, reduction = "harmony", verbose = F) |>
      FindNeighbors(dims = 1:pcs, reduction = "harmony", verbose = F) |>
      FindClusters(resolution = res, reduction = "harmony", verbose = F)
    return(seurat.obj)
  }else{
    seurat.obj <- seurat.obj |>
      RunUMAP(dims = 1:pcs, reduction = "pca") |>
      FindNeighbors(dims = 1:pcs, reduction = "pca") |>
      FindClusters(resolution = res, reduction = "pca")
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
  props$names <- rownames(props) 
  active.idents <- as.data.frame(table(seurat.obj@active.ident, seurat.obj@meta.data[[subject.slot]]))
  
  # exclude sn cell types not used in decon
  pops.decon <- colnames(decon.obj$Var.prop)
  pops.sn <- levels(active.idents$Var1)
  active.idents <- active.idents[active.idents$Var1 %in% pops.decon,]
  active.idents <- active.idents  |>
    group_by(Var2)|>
    mutate(proportion = Freq/sum(Freq))
  count <- nrow(props)
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

# you need to finish the last part of this ratio loop to name the output 