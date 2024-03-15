# Visualize DE 
libs <- c("tidyverse", "wesanderson", "ggrepel") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load the results 
deseq.res <- read.csv("data/processed/models/adjusted_de.csv")
pal <- wes_palette(name = "Zissou1", type = "continuous")
phenotypes <- read.csv("data/processed/bulk/pheno_table.csv")
bulk <- read.csv("data/processed/bulk/all_counts.csv", row.names = 1, check.names = F)
# Set reference factor levels for phenotypes
pheno.reorder <- phenotypes |> 
  mutate(new.id = as.factor(new.id)) |> 
  mutate(Genotype = factor(case_when(
    str_detect(new.id, "WT") ~ "WT",
    str_detect(new.id, "cmAKO") ~ "cmAKO")),
    Treatment = factor(case_when(
      str_detect(new.id, "Sham") ~ "Sham",
      str_detect(new.id, "CAD") ~ "CAD")))

# Unadjusted results first 

# Find top DE genes
top.raw <- deseq.res |> 
  arrange(padj.raw) |>
  slice_head(n = 10)

# Add conditional color formatting for significance
deseq.res <- deseq.res |> 
  mutate(significant = case_when(
    padj.raw >  0.05 ~ FALSE,
    padj.raw <= 0.05 ~ TRUE,
    .default = FALSE
  ))

plot.raw <- ggplot(deseq.res, aes(x = log2FoldChange.raw, y = -log10(padj.raw), color = significant)) +
  geom_point(alpha = 1, size = 6) + 
  scale_color_manual(values = c("#999999", "#ed9209")) +
  geom_text_repel(data = top.raw, aes(label = gene), size = 12, box.padding = unit(0.35, "lines"), 
                  force = 20, segment.linetype = 2, segment.size = 0.6, color = "black", force_pull = 0.01, min.segment.length = 0) + # label most sig genes 
  theme_minimal() +
  labs(x = "log2(Fold Change)", y = "-log10(adjusted p-value)", color = "p < 0.05") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "solid") + # add a line for p=0.05
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black", size = 25),
        axis.ticks = element_blank(),
        title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.key.size = unit(0.7, 'in'),
        legend.title = element_text(size = 28, vjust = 0.7),
        axis.title = element_text(color = "black", size = 28),
        panel.background = element_rect(color="black"),
        plot.background = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"))

# Save 
png(file = "results/10_plot_de/volcano_unadjusted.png",
    width = 12, 
    height = 9,
    units = "in",
    res = 300)

plot.raw

dev.off()

# Repeat for composition adjusted results
# Find top DE genes
top.clr <- deseq.res |> 
  arrange(padj.clr) |>
  slice_head(n = 10)

# Find limits for color scale/legend
limit <- deseq.res |> 
  na.omit() |>
  arrange(desc(abs(p_diff))) |>
  slice_head(n = 1) |>
  pull(p_diff) * c(-1, 1)

plot.clr <- ggplot(deseq.res, aes(x = log2FoldChange.clr, y = -log10(padj.clr), color = p_diff)) +
  geom_point(alpha = 1, size = 6) +
  scale_color_distiller(type = "div", palette = "Spectral" ) +
  geom_text_repel(data = top.clr, aes(label = gene), size = 12, box.padding = unit(0.35, "lines"), 
                  force = 20, segment.linetype = 2, segment.size = 0.6, color = "black", force_pull = 0.01, min.segment.length = 0) + # label most sig genes 
  theme_minimal() +
  labs(x = "log2(Fold Change)", y = "-log10(adjusted p-value)", color = "Change in -log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "solid") + # add a line for p=0.05
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black", size = 25),
        axis.ticks = element_blank(),
        title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.key.size = unit(0.7, 'in'),
        legend.title = element_text(size = 28, vjust = 0.7),
        axis.title = element_text(color = "black", size = 28),
        panel.background = element_rect(color="black"),
        plot.background = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"))

# Save 
png(file = "results/10_plot_de/volcano_adjusted.png",
    width = 12, 
    height = 9,
    units = "in",
    res = 300)

plot.clr

dev.off()



# Plot the PCA of samples
# Run PCA
pca <- prcomp(bulk)

PoV <- pca$sdev^2/sum(pca$sdev^2) * 100  
PoV <- round(PoV, digits = 1)
# Merge with sample info, then plot
pca <- pca$rotation |> 
  as.data.frame() |> 
  dplyr::select(PC1, PC2) 
pca$new.id <- row.names(pca)

my_palette <- c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00")
legend.names <- c("Sham_1","Sham_2", "TAC_1", "TAC_2")

pca <- pca |> left_join(pheno.reorder) |> 
  mutate(gene_treat = paste(genotype, treatment)) |> 
  filter(type == "whole") 
pca.plot <- pca |> 
  ggplot(aes(x = PC1, y = PC2, color = gene_treat)) +
  geom_point(size = 8, color = "black") +
  geom_point(size = 7) +
  
  geom_label_repel(aes(fill = gene_treat),
                   label = pca$new.id, color = "black", alpha = 0.8,
                   box.padding = 0.5, segment.curvature = -0.2,
                   segment.ncp = 3, segment.angle = 20, force = 10, 
                   max.overlaps = 10, force_pull = 0.01, size = 6,
                   nudge_x = 0.015, nudge_y = 0.02) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(expand = expansion(mult = 0.3), name = paste0("PC1", " (", PoV[1], " % of total variance)")) +
  scale_y_continuous(expand = expansion(mult = 0.3), name = paste0("PC2", " (", PoV[2], " % of total variance)")) +
  theme(axis.text.x = element_text(vjust = 0.5),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'),
        plot.margin = unit(c(1,1,1,1), units = "cm"),
        text = element_text(size = 25))
# Save plot to results 
png(file = "results/10_plot_de/pca.png",
    width = 12, 
    height = 9,
    units = "in",
    res = 300)

pca.plot

dev.off()

