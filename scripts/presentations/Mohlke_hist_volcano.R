# Visualize DE 
libs <- c("tidyverse", "wesanderson", "ggrepel", "DESeq2","patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load the results 
# Load the results and expression matrices
res <- readRDS(paste0("data/processed/models/adjusted_de_interaction.RDS")) 
names <- resultsNames(res)[-1]

# make a list of each result
res <- lapply(names, function(x){
  results(res, name=x) |> 
    as.data.frame() |> 
    rownames_to_column(var = "gene")
})
names(res) <- names


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
      str_detect(new.id, "MI") ~ "MI")))
#####
# Patch to plot from list of results contrasts 
####
plotDE <- function(x, title){
  # Find top DE genes
  top <- x |> 
    arrange(padj) |>
    slice_head(n = 10)
  
  # Add conditional color formatting for significance
  x <- x |> 
    mutate(significant = case_when(
      padj >  0.05 ~ FALSE,
      padj <= 0.05 & abs(log2FoldChange) >= 0.583  ~ TRUE,
      .default = FALSE
    ))
  
  plot.fib <- ggplot(x, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 1, size = 6) + 
    scale_color_manual(values = c("#999999", "#ed9209")) +
    geom_text_repel(data = top, aes(label = gene), size = 12, box.padding = unit(0.35, "lines"), 
                    force = 20, segment.linetype = 2, segment.size = 0.6, color = "black", force_pull = 0.01, min.segment.length = 0) + # label most sig genes 
    theme_minimal() +
    xlim(-5, 5) +
    labs(x = "log2(Fold Change)", y = "-log10(adjusted p-value)", 
         color = "p < 0.05 and\nFold Change > 1.5", title = title) +
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
}

p.de <- lapply(1:length(res), function(n){plotDE(res[[n]], names[[n]])})

# Save 
if(!dir.exists("results/10_plot_de")){
  dir.create("results/10_plot_de")
}
png(file = "results/10_plot_de/volcano_clr.png",
    width = 36, 
    height = 18,
    units = "in",
    res = 600)

wrap_plots(p.de)

dev.off()


#####
# Unadjusted results first 

# Find top DE genes
top.raw <- deseq.res |> 
  filter(gene %in% c("Ddr2", "Errfi1", "Egr1", "Pik3r1", "Zbtb16"))

# Add conditional color formatting for significance
deseq.res <- deseq.res |> 
  mutate(significant = case_when(
    padj.raw >  0.05 ~ FALSE,
    padj.raw <= 0.05 & abs(log2FoldChange.raw) >= 0.583  ~ TRUE,
    .default = FALSE
  ))

plot.raw <- ggplot(deseq.res, aes(x = log2FoldChange.raw, y = -log10(padj.raw), color = significant)) +
  geom_point(alpha = 1, size = 6) + 
  scale_color_manual(values = c("#999999", "#ed9209")) +
  geom_text_repel(data = top.raw, aes(label = gene), size = 12, box.padding = unit(0.35, "lines"), 
                  force = 20, segment.linetype = 2, segment.size = 0.6, color = "black", force_pull = 0.01, min.segment.length = 0) + # label most sig genes 
  theme_minimal() +
  xlim(-5, 5) +
  labs(x = "log2(Fold Change)", y = "-log10(adjusted p-value)", color = "p < 0.05 and\nFold Change > 1.5") +
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
if(!dir.exists("results/10_plot_de")){
  dir.create("results/10_plot_de")
}
png(file = "results/mohlke/ihc_volcano_unadj.png",
    width = 12, 
    height = 9,
    units = "in",
    res = 300)

plot.raw

dev.off()

# Repeat for composition adjusted results
# Find top DE genes
top.clr <- deseq.res |> 
  filter(gene %in% c("Ddr2", "Errfi1", "Egr1", "Pik3r1", "Zbtb16")) 

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
  xlim(-6, 6) +
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
png(file = "results/mohlke/ihc_volcano.png",
    width = 12, 
    height = 9,
    units = "in",
    res = 300)

plot.clr

dev.off()
