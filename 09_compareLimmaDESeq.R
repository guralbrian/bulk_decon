# Compare results from DESeq2 and Limma-voom
library(DESeq2)
library(edgeR)
library(tidyverse)
# ? For genes sig in both comparison, do they agree in the direction of effect?
# pre and post adjustment

# load in DESeq2 and Limma-voom results

limma_treat <- read.csv("data/processed/models/limma-treatContrast.csv")
limma_interact <- read.csv("data/processed/models/limma-geneXtreatContrast.csv")

deseq_adj <- readRDS("data/processed/models/adjusted_de_interaction.RDS")
deseq_unadj <- readRDS("data/processed/models/unadjusted_de_interaction.RDS")


#Peyton told me to shrink the DESeq lfc
shrink_inter_adj <- lfcShrink(deseq_adj, coef = 'treatmentMI.genotypecmAKO', type = 'apeglm')
shrink_inter_unadj <- lfcShrink(deseq_unadj, coef = 'treatmentMI.genotypecmAKO', type = 'apeglm')

# Pull interaction effects from DESeq2 results
# Merge and format them to match limma
deseq_join_inter <- merge(
  as.data.frame(shrink_inter_adj), 
  as.data.frame(shrink_inter_unadj), 
  by = "row.names", 
  suffixes = c(".adj", ".unadj")
) |>
  dplyr::rename(gene = Row.names) |>
  mutate(
    p_diff = -log10(padj.adj) - -log10(padj.unadj),
    adj.sig = padj.adj <= 0.05 & padj.unadj >= 0.05,
    unadj.sig = padj.adj >= 0.05 & padj.unadj <= 0.05,
    both.sig = padj.adj <= 0.05 & padj.unadj <= 0.05
  )

# Repeat with Limma results
limma_join_inter <- limma_interact |>
  mutate(
    p_diff = -log10(adj.P.Val.adj) - -log10(adj.P.Val),
    adj.sig = adj.P.Val.adj <= 0.05 & adj.P.Val >= 0.05,
    unadj.sig = adj.P.Val.adj >= 0.05 & adj.P.Val <= 0.05,
    both.sig = adj.P.Val.adj <= 0.05 & adj.P.Val <= 0.05
  )

# Join the two
deseq_simple <- deseq_join_inter |> 
  select(gene,log2FoldChange.unadj, log2FoldChange.adj, padj.unadj, padj.adj)
colnames(deseq_simple) <- c("gene", "lfc-pre", "lfc-post", "p-pre", "p-post")

limma_simple <- limma_join_inter |> 
  select(gene, logFC, logFC.adj, adj.P.Val, adj.P.Val.adj, AveExpr)
colnames(limma_simple) <- c("gene", "lfc-pre", "lfc-post", "p-pre", "p-post", "AveExpr")

interaction_joined <- merge(
  as.data.frame(deseq_simple), 
  as.data.frame(limma_simple), 
  by = "gene", 
  suffixes = c("-de", "-lim")
) 



# Helper functions
sign_agreement <- function(x, y) {
  sum(sign(x) == sign(y)) / length(x)
}

categorize_significance <- function(pvalue) {
  ifelse(pvalue < 0.05, "Significant", "Not Significant")
}

# Compute agreement metrics
metrics <- interaction_joined %>%
  mutate(
    sig_pre_de = categorize_significance(`p-pre-de`),
    sig_post_de = categorize_significance(`p-post-de`),
    sig_pre_lim = categorize_significance(`p-pre-lim`),
    sig_post_lim = categorize_significance(`p-post-lim`),
    
    lfc_agree_pre = sign(`lfc-pre-de`) == sign(`lfc-pre-lim`),
    lfc_agree_post = sign(`lfc-post-de`) == sign(`lfc-post-lim`),
    
    sig_agree_pre = sig_pre_de == sig_pre_lim,
    sig_agree_post = sig_post_de == sig_post_lim
  ) %>%
  summarise(
    lfc_agreement_pre_pct = mean(lfc_agree_pre, na.rm = F) * 100,
    lfc_agreement_post_pct = mean(lfc_agree_post, na.rm = F) * 100,
    
    sig_agreement_pre_pct = mean(sig_agree_pre, na.rm = T) * 100,
    sig_agreement_post_pct = mean(sig_agree_post, na.rm = T) * 100,
    
    sig_change_agreement = mean((sig_pre_de != sig_post_de) == (sig_pre_lim != sig_post_lim), na.rm = T) * 100
  )


metrics
# Tally significance 
# Function to categorize significance

# Add significance categories for each test and adjustment
interaction_joined <- interaction_joined %>%
  mutate(
    sig_pre_de = categorize_significance(`p-pre-de`),
    sig_post_de = categorize_significance(`p-post-de`),
    sig_pre_lim = categorize_significance(`p-pre-lim`),
    sig_post_lim = categorize_significance(`p-post-lim`)
  )

# Count genes with changed significance pre vs post adjustment for DESeq2 and limma
sig_changes <- interaction_joined %>%
  summarise(
    DESeq2_sig_loss = sum(sig_pre_de == "Significant" & sig_post_de == "Not Significant", na.rm = T),
    DESeq2_sig_gain = sum(sig_pre_de == "Not Significant" & sig_post_de == "Significant", na.rm = T),
    limma_sig_loss = sum(sig_pre_lim == "Significant" & sig_post_lim == "Not Significant", na.rm = T),
    limma_sig_gain = sum(sig_pre_lim == "Not Significant" & sig_post_lim == "Significant", na.rm = T)
  )
print(sig_changes)


# Prepare data correctly for scatter plot
plot_data <- interaction_joined %>%
  transmute(
    gene,
    DESeq2_pre = -log10(`p-pre-de`),
    DESeq2_post = -log10(`p-post-de`),
    Limma_pre = -log10(`p-pre-lim`),
    Limma_post = -log10(`p-post-lim`)
  ) %>%
  pivot_longer(cols = -gene, names_to = c("Method", "Adjustment"), names_sep = "_", values_to = "neg_log10_p") %>%
  pivot_wider(names_from = Method, values_from = neg_log10_p)

# Scatter plot showing pre and post adjustments in different colors
ggplot(plot_data, aes(x = DESeq2, y = Limma, color = Adjustment)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_abline(slope = 1, color = "black") +
  labs(x = "-log10 p-value DESeq2", y = "-log10 p-value Limma", color = "Adjustment",
       title = "Comparison of p-values for DESeq2 vs Limma") +
  theme_minimal() +
  #xlim(c(0,5)) +
  #ylim(c(0,5)) +
  scale_color_manual(values = c("pre" = "steelblue", "post" = "darkorange"))


# Prepare lfc data 
plot_lfc <- interaction_joined %>%
  transmute(
    gene,
    DESeq2_pre = `lfc-pre-de`,
    DESeq2_post = `lfc-post-de`,
    Limma_pre = `lfc-pre-lim`,
    Limma_post = `lfc-post-lim`
  ) %>%
  pivot_longer(cols = -gene, names_to = c("Method", "Adjustment"), names_sep = "_", values_to = "lfc") %>%
  pivot_wider(names_from = Method, values_from = lfc)

# Scatter plot showing pre and post adjustments in different colors
ggplot(plot_lfc, aes(x = DESeq2, y = Limma, color = Adjustment)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(slope = 1, color = "black") +
  labs(x = "logFC DESeq2", y = "logFC Limma", color = "Adjustment",
       title = "Comparison of logFC for DESeq2 vs Limma") +
  theme_minimal() +
  scale_color_manual(values = c("pre" = "steelblue", "post" = "darkorange"))


# Look at different parts of the objects used in the DESeq and limma scripts
# These lines only work if you have the intermediatary objecst from each script in your environment
joined_norm_factors <- data.frame(sample = names(sizeFactors(dds.clr)),
                                  deseq_sizeFactors = sizeFactors(dds.clr),
                                  limma_normFactors = d0$samples$norm.factors)
joined_norm_factors |> 
  ggplot(aes(x = deseq_sizeFactors, y = limma_normFactors)) +
  geom_text(aes(label = sample)) +
  geom_abline(slope = 1) +
  theme_minimal()

cor(sample.clr$clr.Cardiomyocytes, sizeFactors(dds.clr))
cor(sample.clr$clr.Fibroblast, sizeFactors(dds.clr))
cor(sample.clr$clr.Cardiomyocytes, d0$samples$norm.factors)
cor(sample.clr$clr.Fibroblast, d0$samples$norm.factors)
plot(sample.clr$clr.Cardiomyocytes, counts["Fhl1", ], 
     xlab="CLR Cardiomyocytes", ylab="Fhl1 Expression")
