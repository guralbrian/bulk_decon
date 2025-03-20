# Limma Voom version of DESeq2 pipline

library(edgeR)
library(tidyverse)

# Load and format input data
decon.whole <- read.csv("data/processed/compositions/whole_samples.csv")
phenotypes <- read.csv("data/processed/bulk/pheno_table.csv") |> 
  filter(type == "whole")
bulk.all <- read.csv("data/processed/bulk/all_counts.csv", row.names = 1, check.names = FALSE)

# Format bulk data and compositions
bulk.all <- bulk.all[,unique(decon.whole$new.id)]
counts <- mutate_all(bulk.all, \(x) round(as.numeric(as.character(x))))

# Reorder the data frame
counts <- counts[,phenotypes$new.id]

# filter out genes with less than 10 counts in 75% of samples
keep <- rowSums(counts >= 10) >= 4

# Make edgeR DGEList object and calc normalization factors (accounts for library size)
d0 <- DGEList(counts) |> calcNormFactors()

drop <- which(apply(cpm(d0), 1, max) < 1)
d <- d0[-drop,] 
dim(d)
# Make group variables

# Prepare sample information
sample_info <- phenotypes |>
  mutate(across(c( genotype, treatment), as.factor))

sample_info$treatment <- relevel(sample_info$treatment, ref = "Sham")
sample_info$genotype <- relevel(sample_info$genotype, ref = "WT")
sample.clr <- cbind(sample_info, comps.clr[colnames(bulk),])
# define the vars as their own objects
genotype <- relevel(sample_info$genotype, ref = "WT")
treatment <- relevel(sample_info$treatment, ref = "Sham")


#### Voom transformation and calc of variance weights
mm <- model.matrix(~genotype * treatment)

y <- voom(d, mm, plot = T)

# Fit a model for each gene using weighted least squares
fit <- lmFit(y, mm)

# Comparisons/logFC are obtained as contrasts of the fitted lm
tmp <- contrasts.fit(fit, coef = 4)

#shrinks std err that are much far from the average std err  
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf)

# Repeat with the cell-type covariates
library(compositions)

# Prepare CLR transformed compositions
decon.wide <- decon.whole |> 
  dplyr::select(new.id, CellType, Prop) |> 
  pivot_wider(names_from = "new.id", values_from = "Prop") |> 
  column_to_rownames("CellType") |>
  mutate(across(everything(), as.numeric))

comps <- compositions::clr(t(decon.wide)) |> as.data.frame()
colnames(comps) <- paste0("clr.", colnames(comps)) 

# Get vectors of compositions
cardiomyocytes <- comps.clr$clr.Cardiomyocytes
fibroblasts <- comps.clr$clr.Fibroblast


# Voom transformation and calc of variance weights
mm <- model.matrix(~genotype * treatment + cardiomyocytes + fibroblasts)

y <- voom(d, mm, plot = T)

# Fit a model for each gene using weighted least squares
fit.adj <- lmFit(y, mm)

# Comparisons/logFC are obtained as contrasts of the fitted lm
tmp <- contrasts.fit(fit.adj, coef = 6)

#shrinks std err that are much far from the average std err  
tmp <- eBayes(tmp)

top.table.adj <- topTable(tmp, sort.by = "P", n = Inf)
colnames(top.table.adj) <- paste0(colnames(top.table.adj), ".adj") 
top.table.adj$gene <- row.names(top.table.adj) 

# Join the results for the interaction contrast
top.table$gene <- row.names(top.table)

limma.inter.res <- left_join(top.table, top.table.adj)

write.csv(limma.inter.res, "data/processed/models/limma-geneXtreatContrast.csv", row.names = FALSE)

##### Repeat everything for the treatment contrast
# Unadjusted model
tmp <- contrasts.fit(fit, coef = 3)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table$gene <- row.names(top.table)
# Adjusted model
tmp <- contrasts.fit(fit.adj, coef = 3)
tmp <- eBayes(tmp)
top.table.adj <- topTable(tmp, sort.by = "P", n = Inf)
colnames(top.table.adj) <- paste0(colnames(top.table.adj), ".adj") 
top.table.adj$gene <- row.names(top.table.adj) 


limma.inter.res <- left_join(top.table, top.table.adj)

write.csv(limma.inter.res, "data/processed/models/limma-treatContrast.csv", row.names = FALSE)
