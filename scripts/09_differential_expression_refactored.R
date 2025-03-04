# Load libraries
libs <- c("tidyverse", "compositions", "DESeq2")
invisible(lapply(libs, require, character.only = TRUE))

# Load and format input data
decon.whole <- read.csv("data/processed/compositions/whole_samples.csv")
phenotypes <- read.csv("data/processed/bulk/pheno_table.csv") |> 
  filter(type == "whole")
bulk.all <- read.csv("data/processed/bulk/all_counts.csv", row.names = 1, check.names = FALSE)

# Format bulk data and compositions
bulk.all <- bulk.all[,unique(decon.whole$new.id)]
bulk <- mutate_all(bulk.all, \(x) round(as.numeric(as.character(x))))

# Prepare CLR transformed compositions
decon.wide <- decon.whole |> 
  dplyr::select(new.id, CellType, Prop) |> 
  pivot_wider(names_from = "new.id", values_from = "Prop") |> 
  column_to_rownames("CellType") |>
  mutate(across(everything(), as.numeric))

comps.clr <- compositions::clr(t(decon.wide))
colnames(comps.clr) <- paste0("clr.", colnames(comps.clr))

# Prepare sample information
sample_info <- phenotypes |>
  mutate(across(c(new.id, genotype, treatment), as.factor)) |>
  column_to_rownames("new.id")

sample_info$treatment <- relevel(sample_info$treatment, ref = "Sham")
sample_info$genotype <- relevel(sample_info$genotype, ref = "WT")
sample_info <- sample_info[colnames(bulk),]
sample.clr <- cbind(sample_info, comps.clr[colnames(bulk),])

run_deseq <- function(bulk, sample_info, design, output_path) {
  dds <- DESeqDataSetFromMatrix(
    countData = bulk,
    colData = sample_info,
    design = design
  )
  keep <- rowSums(counts(dds) >= 10) >= 4
  dds <- dds[keep,]
  dds <- DESeq(dds)
  saveRDS(dds, output_path)
  dds
}

# Run DESeq2 with and without compositions
dds.clr <- run_deseq(
  bulk, 
  sample.clr,
  ~ treatment + genotype + treatment:genotype + clr.Cardiomyocytes + clr.Fibroblast,
  "data/processed/models/adjusted_de_interaction.RDS"
)

dds.raw <- run_deseq(
  bulk,
  sample_info,
  ~ treatment + genotype + treatment:genotype,
  "data/processed/models/unadjusted_de_interaction.RDS"
)

# Compare results
res.clr <- results(dds.clr, name="treatmentMI.genotypecmAKO")
res.raw <- results(dds.raw, name="treatmentMI.genotypecmAKO")

comparison <- merge(
  as.data.frame(res.clr), 
  as.data.frame(res.raw), 
  by = "row.names", 
  suffixes = c(".clr", ".raw")
) |>
  dplyr::rename(gene = Row.names) |>
  mutate(
    p_diff = -log10(padj.clr) - -log10(padj.raw),
    clr.sig = padj.clr <= 0.05 & padj.raw >= 0.05,
    raw.sig = padj.clr >= 0.05 & padj.raw <= 0.05,
    both.sig = padj.clr <= 0.05 & padj.raw <= 0.05
  )

write.csv(comparison, "data/processed/models/adjusted_de.csv", row.names = FALSE)
