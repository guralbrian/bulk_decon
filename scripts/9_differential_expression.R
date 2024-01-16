# Test for differential gene expression with DESeq2 
# Include composition covariates to test if they mediate DE

# Load libs
libs <- c("tidyverse", "compositions", "reshape2", "stats", "DESeq2", "wesanderson", "ggrepel") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#### Loading and formatting of data ####
# Load compositions, expression, and phenotype data
decon.whole <- read.csv("data/processed/compositions/whole_samples.csv")
phenotypes.real <- read.csv("data/processed/bulk/jensen_pheno.csv")
jensen.bulk <- read.csv("data/processed/bulk/jensen_bulk_clean.csv",  row.names = 1, check.names = F)


# Unmelt for clr transformation 
decon.wide <- decon.whole  |> 
  select(Sub, CellType, Prop) |> 
  pivot_wider(names_from = "Sub", values_from = "Prop") |> 
  column_to_rownames("CellType") %>%
  mutate_all(as.numeric)

# use clr or ilr to incorp the compositions into the design matrix
comps.clr <- compositions::clr(t(decon.wide))

colnames(comps.clr) <- paste0("clr.", colnames(comps.clr))

# Prepare for DESeq2
bulk <- mutate_all(jensen.bulk, function(x) round(as.numeric(as.character(x)), digits = 0)) # round to integers

# Set reference factor levels for phenotypes
pheno.reorder <- phenotypes.real |> 
  mutate(Sub = as.factor(de_id)) |> 
  mutate(Genotype = factor(case_when(
    str_detect(Sub, "WT") ~ "WT",
    str_detect(Sub, "KO") ~ "KO")),
    Treatment = factor(case_when(
      str_detect(Sub, "sham") ~ "Sham",
      str_detect(Sub, "MI") ~ "TAC")))

# Prepare sample information
sample_info <- data.frame(
  row.names = colnames(bulk),
  genotype = as.factor(pheno.reorder$Genotype),
  treatment = as.factor(pheno.reorder$Treatment)
)

sample_info$treatment <- relevel(sample_info$treatment, ref = "Sham")
sample_info$genotype <- relevel(sample_info$genotype, ref = "KO")

# Add clr transforms to sample info
sample.clr <- cbind(sample_info, comps.clr)

#### Run DESeq2 ####
## DESeq with compositions
# Create a DESeqDataSet
dds.clr <- DESeqDataSetFromMatrix(
  countData = bulk,
  colData = sample.clr,
  design = ~ treatment + genotype + treatment:genotype + clr.Cardiomyocytes + clr.Fibroblast
)

# Run DESeq
dds.clr <- DESeq(dds.clr)

# Pull out interaction term results
res.clr <- results(dds.clr, name="treatmentTAC.genotypeWT")

## DESeq without compositions
# Create a DESeqDataSet
dds.raw <- DESeqDataSetFromMatrix(
  countData = bulk,
  colData = sample_info,
  design = ~ treatment + genotype + treatment:genotype
)

# Run DESeq 
dds.raw <- DESeq(dds.raw)

# Pull out interaction term results
res.raw <- results(dds.raw, name="treatmentTAC.genotypeWT")

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