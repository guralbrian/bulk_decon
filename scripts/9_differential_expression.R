# Test for differential gene expression with DESeq2 
# Include composition covariates to test if they mediate DE

# Load libs
libs <- c("tidyverse", "compositions", "reshape2", "stats", "DESeq2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#### Loading and formatting of data ####
# Load compositions, expression, and phenotype data
decon.whole <- read.csv("data/processed/compositions/whole_samples.csv")
phenotypes <- read.csv("data/processed/bulk/pheno_table.csv")
bulk.all <- read.csv("data/processed/bulk/all_counts.csv", row.names = 1, check.names = F)

# Subset to whole samples
bulk.all <- bulk.all[,c(unique(decon.whole$new.id)),]
phenotypes <- phenotypes |> 
  filter(type == "whole")

# Unmelt for clr transformation 
decon.wide <- decon.whole  |> 
  dplyr::select(new.id, CellType, Prop) |> 
  pivot_wider(names_from = "new.id", values_from = "Prop") |> 
  column_to_rownames("CellType") %>%
  mutate_all(as.numeric)
#decon.wide[decon.wide==0] <- 0.01

# use clr or ilr to incorp the compositions into the design matrix
comps.clr <- compositions::clr(t(decon.wide))

colnames(comps.clr) <- paste0("clr.", colnames(comps.clr))

# Prepare for DESeq2
bulk <- mutate_all(bulk.all, function(x) round(as.numeric(as.character(x)), digits = 0)) # round to integers

# Set reference factor levels for phenotypes
pheno.reorder <- phenotypes |> 
  mutate(new.id = as.factor(new.id),
         genotype = as.factor(genotype),
         treatment = as.factor(treatment))

# Prepare sample information
sample_info <- data.frame(
  row.names = pheno.reorder$new.id,
  genotype = pheno.reorder$genotype,
  treatment = pheno.reorder$treatment
)

sample_info$treatment <- relevel(sample_info$treatment, ref = "Sham")
sample_info$genotype <- relevel(sample_info$genotype, ref = "WT")
sample_info <- sample_info[colnames(bulk),]

# Add clr transforms to sample info
sample.clr <- cbind(sample_info, comps.clr[colnames(bulk),])

#### Run DESeq2 ####
## DESeq with compositions
# Create a DESeqDataSet
dds.clr <- DESeqDataSetFromMatrix(
  countData = bulk,
  colData = sample.clr,
  design = ~ treatment + genotype + treatment:genotype + clr.Fibroblast + clr.Cardiomyocytes #+ clr.Macrophage
)

# Filter out lowly expressed genes
smallestGroupSize <- 4
keep <- rowSums(counts(dds.clr) >= 10) >= smallestGroupSize
dds.clr <- dds.clr[keep,]

# remove outlier
keep <- colnames(dds.clr)[colnames(dds.clr) != "WT Sham (4)"]
dds.clr <- dds.clr[,keep]

# Run DESeq
dds.clr <- DESeq(dds.clr)

# Save the results 
#Save outputs
if(!dir.exists("data/processed/models")){
  dir.create("data/processed/models")
}
saveRDS(dds.clr, "data/processed/models/adjusted_de_interaction.RDS")

# Pull out interaction term results
res.clr <- results(dds.clr, name="treatmentMI.genotypecmAKO")

## DESeq without compositions
# Create a DESeqDataSet
dds.raw <- DESeqDataSetFromMatrix(
  countData = bulk,
  colData = sample_info,
  design = ~ treatment + genotype + treatment:genotype
)

# Filter out lowly expressed genes
smallestGroupSize <- 4
keep <- rowSums(counts(dds.raw) >= 15) >= smallestGroupSize
dds.raw <- dds.raw[keep,]

# remove outlier
keep <- colnames(dds.clr)[colnames(dds.clr) != "WT Sham (4)"]
dds.clr <- dds.clr[,keep]

# Run DESeq 
dds.raw <- DESeq(dds.raw)

# Save the results 
saveRDS(dds.raw, "data/processed/models/unadjusted_de_interaction.RDS")

# Pull out interaction term results
res.raw <- results(dds.raw, name="treatmentMI.genotypecmAKO")

#### Contrast DE results between clr and raw/unadjusted runs ####
# Compare results
comparison.clr <- merge(as.data.frame(res.clr), as.data.frame(res.raw), by = "row.names", suffixes = c(".clr", ".raw"))
colnames(comparison.clr)[1] <- "gene"
# Calculate difference in -log10 p-values
comparison.clr$p_diff <- -log10(comparison.clr$padj.clr) -  # low numbers mean that its not signficant when clr is used
  -log10(comparison.clr$padj.raw)    # large numbers mean that its significant w/o comp

# Make categorical column
comparison.clr <- comparison.clr |>
  mutate(clr.sig = ifelse(padj.clr <= 0.05 & padj.raw >= 0.05, T, F),
         raw.sig = ifelse(padj.clr >= 0.05 & padj.raw <= 0.05, T, F),
         both.sig = ifelse(padj.clr <= 0.05 & padj.raw <= 0.05, T, F)) 

# Save the results 
write.csv(comparison.clr, "data/processed/models/adjusted_de.csv", row.names = F)
