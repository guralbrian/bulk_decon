# sn -> expression profiles for each cell type 
# dataframe to organize cell type ratios for each sample
# matrix of expression of each sample
# script to load the matrix and run DESeq on a subset of samples

# Script to be paired with shell script to run many versions in parallel

# Load Libraries

# List libraries
libs <- c("compositions", "DESeq2", "stringr", "dplyr") # list libraries here
# Require all of them
lapply(libs, require, character.only = T)

rm(libs)

# Get commandArgs
args <- commandArgs(trailingOnly = TRUE)
model_arg <-  args[[1]] |> as.numeric()
gene_batch <- args[[2]] |> as.numeric()


# Gene batch size
batch_size <- 500

# Calculate start and end row for the chunk
start_row <- (gene_batch - 1) * batch_size + 1
end_row <- gene_batch * batch_size

# Read and subset the ratios and counts matrix
sim.counts <- read.csv("data/processed/deseq_simulation/simulated_counts.csv", row.names = 1)[c(start_row:end_row),]
#genes.downsampled <- sample(row.names(sim.counts), 2000, replace = F)
#sim.counts <- sim.counts[genes.downsampled,]
ratios <- read.csv("data/processed/deseq_simulation/simulated_ratios.csv", row.names = 1)
rownames(ratios) <- ratios$sample
cell.types <- colnames(ratios)[!(colnames(ratios) %in% c("pct.change", "sample"))]

# Make list of models
models <- list(
  unadjusted    = ~ 0 + pct.change,
  raw_cm        = ~ 0 + pct.change + Cardiomyocytes,
  raw_cm_fb     = ~ 0 + pct.change + Cardiomyocytes + Fibroblast,
  clr_cm        = ~ 0 + pct.change + clr_Cardiomyocytes,
  clr_cm_fb     = ~ 0 + pct.change + clr_Cardiomyocytes + clr_Fibroblast,
  pc1           = ~ 0 + pct.change + PC1,
  pc2           = ~ 0 + pct.change + PC1 + PC2
  )

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

sample_info$pct.change <- sample_info$pct.change |>
  str_replace_all("-", "_") |>
  as.factor()

# Relevel the sample info to have the no change group as the reference
sample_info$pct.change <- relevel(sample_info$pct.change, ref = "0")

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

# Run the DESeq 
dds <- DESeq(dds, fitType="mean")

# Get the differential expression analysis results, contrasted against the 0 percent change group
genes <- lapply(pct.use, function(x){
  results(dds, contrast = c("pct.change", "0", x)) |> as.data.frame()})

# rename the nested list with the percent changes
names(genes) <- pct.use #lapply(str_split(resultsNames(dds)[-1], "pct.change"), "[[", 2) |> unlist()


sig.genes <- lapply(genes, function(x){
  sig.genes.temp <- x |> 
    mutate(sig = case_when(
      padj <= 0.05 ~ TRUE,
      padj > 0.05 ~ FALSE
    )) |> 
    pull(sig) |> 
    table()
  if(length(sig.genes.temp) == 1){
    sig.genes.temp[[2]] <- 0
  }
  pct.sig <- sig.genes.temp[[2]]/sum(sig.genes.temp)
  return(pct.sig)
})

sig.df <- data.frame(pct.change = names(sig.genes),
                     pct.sig = unlist(sig.genes),
                     model = model_arg,
                     batch = gene_batch)

if(!dir.exists("data/processed/deseq_simulation/batched_output/")){
  dir.create("data/processed/deseq_simulation/batched_output/")
}

write.csv(sig.df, 
          paste0("data/processed/deseq_simulation/batched_output/",model_arg, "_",gene_batch, ".csv"), 
          row.names = F)
