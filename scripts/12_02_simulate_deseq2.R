# sn -> expression profiles for each cell type 
# dataframe to organize cell type ratios for each sample
# matrix of expression of each sample
# script to load the matrix and run DESeq on a subset of samples

# Script to be paired with shell script to run many versions in parallel

# Load Libraries

# List libraries
libs <- c("compositions", "DESeq2", "stringr", "tidyverse") # list libraries here
# Require all of them
lapply(libs, require, character.only = T)

rm(libs)

# Get commandArgs
args <- commandArgs(trailingOnly = TRUE)
model_arg <-  args[[1]] |> as.numeric()
gene_batch <- args[[2]] |> as.numeric()


# Gene batch size
batch_size <- 1000

# Calculate start and end row for the chunk
start_row <- (gene_batch - 1) * batch_size + 1
end_row <- gene_batch * batch_size

# Load list of true positives
counts.dict <-  read_lines("data/processed/deseq_simulation/true_postives.txt")

# Read and subset the ratios and counts matrix
sim.counts <- read.csv("data/processed/deseq_simulation/simulated_counts.csv", row.names = 1)[c(start_row:end_row),]

ratios <- read.csv("data/processed/deseq_simulation/simulated_ratios.csv", row.names = 1)
rownames(ratios) <- ratios$sample
cell.types <- colnames(ratios)[!(colnames(ratios) %in% c("pct.change", "sample"))]

# Make list of models
models <- list(
  unadjusted    = ~ 0 + pct.change,
  raw_cm        = ~ 0 + pct.change + Cardiomyocytes,
  raw_cm_fb     = ~ 0 + pct.change + Cardiomyocytes + Fibroblast,
  clr_cm        = ~ 0 + pct.change + clr_Cardiomyocytes,
  clr_cm_fb     = ~ 0 + pct.change +  clr_Cardiomyocytes + clr_Fibroblast,
  pc1           = ~ 0 + pct.change + PC1
  )

var_names <- c("Cardiomyocytes", "Fibroblast", "clr_Cardiomyocytes", "clr_Fibroblast", "PC1")
# Compute PCA and clr 
clr.cells <- clr(ratios[,cell.types])
colnames(clr.cells) <- paste0("clr_", colnames(clr.cells))
pca.cells <- prcomp(ratios[,cell.types])$x

# Define what model we'll use
model.use <- models[[model_arg]]
# Change formatting so that DESeq doesn't complain
sample_info <- ratios[colnames(sim.counts),] 
sample_info <- merge(sample_info, clr.cells, by = "row.names")
sample_info <- merge(sample_info, pca.cells, by.x = "Row.names", by.y = "row.names")

# There's an issue converting the numbers to characters, so we need to us this weird way
# .. its bc going to character is converting the whole binary representation
sample_info$pct.change <- sprintf("%.2f", sample_info$pct.change) |> 
  str_replace_all("-", "_") |> 
  as.factor()
# Relevel the sample info to have the no change group as the reference
sample_info$pct.change <- relevel(sample_info$pct.change, ref = "0.00")

row.names(sample_info) <- sample_info$Row.names
sample_info <- sample_info |> select(-Row.names)
sample_info <- sample_info[colnames(sim.counts),]

# Get the list of group percents to use
pct.use <- levels(sample_info$pct.change)[-1]

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = sim.counts,
  colData = sample_info,
  design = model.use
)

# Run DESeq 
dds <- DESeq(dds)


# Get the differential expression analysis results, contrasted against the 0 percent change group
genes <- lapply(pct.use, function(x){
  results(dds, contrast = c("pct.change", "0.00", x)) |> as.data.frame()})

# rename the nested list with the percent changes
names(genes) <- pct.use 

non.pct.vars <- resultsNames(dds)[!str_detect(resultsNames(dds), "pct.change")]

# add the cell type covariate if it exists
for(i in non.pct.vars){
  genes[[i]] <- results(dds, name = i) |> as.data.frame()
}

sig.genes <- lapply(1:length(genes), function(x){
  genes.temp <- genes[[x]]
  sig.genes.temp <- genes.temp |> 
    rownames_to_column(var = "gene") |> 
  # Label each gene result by true/false postive/negative
    mutate(
      sig = case_when(  
        padj <= 0.05 ~ TRUE,
        padj > 0.05 ~ FALSE,
        is.na(padj) ~ FALSE),
      true_de = case_when(
        gene %in% counts.dict ~ TRUE,
        .default = FALSE),
      result_category = case_when(
        sig == F & true_de == F ~ "True negative",
        sig == T & true_de == F ~ "False positive",
        sig == F & true_de == T ~ "False negative",
        sig == T & true_de == T ~ "True positive"),
  # Relevel so that all options appear in the resulting list
      result_category = factor(result_category, levels = c("True negative","False positive", 
                                                           "False negative","True positive"))) |> 
    pull(result_category) |> 
    table()
  # make it into a df
  df.temp <- matrix(sig.genes.temp, nrow = 1) |> as.data.frame()
  colnames(df.temp) <- names(sig.genes.temp)
  df.temp$batch <- gene_batch
  df.temp$model <- model_arg
  df.temp$pct.change <- names(genes)[[x]]
  df.temp$nDEGs <- sum(rownames(genes.temp) %in% counts.dict)
  return(df.temp)
})

# merge all results
sig.df <- do.call("rbind", sig.genes)

if(!dir.exists("data/processed/deseq_simulation/batched_output/")){
  dir.create("data/processed/deseq_simulation/batched_output/")
}

write.csv(sig.df, 
          paste0("data/processed/deseq_simulation/batched_output/", model_arg, "_",gene_batch, ".csv"), 
          row.names = F)
