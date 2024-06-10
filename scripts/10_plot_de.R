# Visualize DE 
libs <- c("tidyverse", "ggrepel", "DESeq2","patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load the results and expression matrices
res.adj <- readRDS("data/processed/models/adjusted_de_interaction.RDS")
names <- resultsNames(res.adj)[-1]

# make a list of each result
res.adj <- lapply(names, function(x){
  results(res.adj, name=x) |> 
    as.data.frame() |> 
    rownames_to_column(var = "gene") |> 
    mutate(model = "adjusted")
})
names(res.adj) <- paste0(names, ".adj")

# Repeat for unadjusted results, then combine them
res.unadj <- readRDS("data/processed/models/unadjusted_de_interaction.RDS")
names <- resultsNames(res.unadj)[-1]

# make a list of each result
res.unadj <- lapply(names, function(x){
  results(res.unadj, name=x) |> 
    as.data.frame() |> 
    rownames_to_column(var = "gene")|> 
    mutate(model = "unadjusted")
})

names(res.unadj) <- paste0(names, ".unadj")

# Combine
res <- c(res.adj, res.unadj)

rm(res.unadj)
rm(res.adj)

# make supplement volcanos for each variable/model

dictionary <- data.frame(name = names(res)) |> 
  mutate(
    model = case_when(
      str_detect(name, ".unadj") ~ "unadjusted",
      str_detect(name, ".adj") ~ "adjusted"),
    variable = case_when(
      str_detect(name, "MI_vs_Sham") ~ "MI",
      str_detect(name, "cmAKO_vs_WT") ~ "cmAKO",
      str_detect(name, "Cardiomyocytes") ~ "Cardiomyocytes",
      str_detect(name, "Fibroblast") ~ "Fibroblast",
      str_detect(name, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI")
  )

#####
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
  
  # Make a title for the plot 
  joined.title <- paste0(dictionary[which(dictionary$name == title), "variable"], ", ",
                         dictionary[which(dictionary$name == title), "model"])
  
  plot.fib <- ggplot(x, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 1, size = 6) + 
    scale_color_manual(values = c("#999999", "#ed9209")) +
    geom_text_repel(data = top, aes(label = gene), size = 12, box.padding = unit(0.35, "lines"), 
                    force = 20, segment.linetype = 2, segment.size = 0.6, color = "black", force_pull = 0.01, min.segment.length = 0) + # label most sig genes 
    theme_minimal() +
    xlim(-5, 5) +
    labs(x = "log2(Fold Change)", y = "-log10(adjusted p-value)", 
         color = "p < 0.05 and\nFold Change > 1.5", title = joined.title) +
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

unadj <- res[c(which(dictionary$model == "unadjusted"))]
p.de.unadj <- lapply(1:length(unadj), function(n){plotDE(unadj[[n]], names(unadj)[[n]])})

# Save 
png(file = "results/supp_figs/volcano_unadjusted.png",
    width = 10, 
    height = 30,
    units = "in",
    res = 300)

wrap_plots(p.de.unadj, ncol = 1) 

dev.off()


# Repeat for adjusted results

adj <- res[c(which(dictionary$model == "adjusted"))]
p.de.adj <- lapply(1:length(adj), function(n){plotDE(adj[[n]], names(adj)[[n]])})

# Save 
png(file = "results/supp_figs/volcano_adjusted.png",
    width = 20, 
    height = 30,
    units = "in",
    res = 300)

wrap_plots(p.de.adj, ncol = 2) 

dev.off()


#### Plot unadjusted volcano of interaction term ####

res.unadj <- res$treatmentMI.genotypecmAKO.unadj

# Find top DE genes
top <- res.unadj |> 
  arrange(padj) |>
  slice_head(n = 10)

# Add conditional color formatting for significance
res.unadj <- res.unadj |> 
  mutate(significant = case_when(
    padj >  0.05 ~ FALSE,
    padj <= 0.05 & abs(log2FoldChange) >= 0.583  ~ TRUE,
    .default = FALSE
  ))

p.unadj <- ggplot(res.unadj, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 1, size = 6) + 
  scale_color_manual(values = c("#999999", "#ed9209")) +
  geom_text_repel(data = top, aes(label = gene), size = 12, box.padding = unit(0.35, "lines"), 
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
png(file = "results/10_plot_de/volcano_unadjusted.png",
    width = 12, 
    height = 9,
    units = "in",
    res = 600)

p.unadj

dev.off()

#### Make plot that shows before and after adj ####

res.compare <- full_join(res$treatmentMI.genotypecmAKO.unadj, res$treatmentMI.genotypecmAKO.adj)


# Define your gene list
gene_list <- c("Zbtb16", "Pik3r1", "Egr1")

# Modify the dataset with new color and alpha values
res.compare <- res.compare %>% 
  mutate(
    in_list = gene %in% gene_list,
    color.cat = case_when(
      in_list & model == "unadjusted" ~ paste(gene, "pre-adjustment"),
      in_list & model == "adjusted" ~ paste(gene, "post-adjustment"),
      !in_list & model == "unadjusted" ~ "Unadjusted",
      !in_list & model == "adjusted" ~ "Adjusted",
      TRUE ~ "Other"
    ),
    alpha.val = ifelse(in_list, 1, 0.3)
  )

# Define colors for plotting
colors <- c(
  "unadjusted" = "#fee090", 
  "adjusted" = "#abd9e9" 
)

# Define the subset of data for labeling
res.labeled <- res.compare %>% 
  filter(in_list) 

x.start <- res.labeled |> 
  filter(model == "unadjusted") |> 
  pull(log2FoldChange)
x.end <- res.labeled |> 
  filter(model == "adjusted") |> 
  pull(log2FoldChange)
y.start <- res.labeled |> 
  filter(model == "unadjusted") |> 
  pull(padj)
y.end <- res.labeled |> 
  filter(model == "adjusted") |> 
  pull(padj)

res.labels <- res.labeled %>% 
  filter(model == "adjusted") 

# Create the plot
p.adjusted.de <- ggplot(res.compare, aes(x = log2FoldChange, y = -log10(padj), color = model)) +
  geom_point(alpha = 0.3, size = 5) +  # Base layer for all points with reduced visibility
  geom_point(data = res.labeled, color = "black", size = 10, alpha = 1) +  # Gives colored points a black outline
  geom_point(data = res.labeled, size = 9, alpha = 1) +  # Colored points
  geom_segment(aes(x = x.start[[2]], y = -log10(y.start[[2]])- 0.5, xend = x.end[[2]], yend = -log10(y.end[[2]])+ 0.5), 
               color = "black", arrow = arrow(length = unit(0.5, "cm")), linetype=2) +
  geom_segment(aes(x = x.start[[3]], y = -log10(y.start[[3]])- 0.5, xend = x.end[[3]], yend = -log10(y.end[[3]])+ 0.5), 
               color = "black", arrow = arrow(length = unit(0.5, "cm")), linetype=2) +
  #geom_segment(aes(x = x.start, y = -log10(y.start), xend = x.end, yend = -log10(y.end)), arrow = arrow(length = unit(0.5, "cm"))) +
  geom_text_repel(data = res.labels, aes(label = gene), size = 10, box.padding = unit(0.6, "lines"), 
                  force = 20, color = "black", nudge_x = 0.7, min.segment.length = 5) +  # Labels only for colored points
  scale_color_manual(values = colors, 
                     guide = guide_legend(override.aes = list(alpha = 1, size = 10)),
                     breaks = c("adjusted", "unadjusted")) + 
  theme_minimal() +
  xlim(-5, 5) +
  labs(x = "log2(Fold Change)", y = "-log10(adjusted p-value)", color = "Model Type") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "solid", alpha = 0.5) +
  theme(legend.position = "bottom",
        text = element_text(color = "black", size = 27),
        legend.text = element_text(color = "black", size = 27),
        panel.background = element_rect(color="black"),
        plot.background = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm")) 

# Save 
if(!dir.exists("results/10_plot_de")){
  dir.create("results/10_plot_de")
}
png(file = "results/10_plot_de/volcano_adjusted.png",
    width = 12, 
    height = 9,
    units = "in",
    res = 600)

p.adjusted.de

dev.off()






# Plot the PCA of samples
# Run PCA
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

# Convert to CPM
cpm <- apply(bulk,2, function(x) (x/sum(x))*1000000) |> 
  as.data.frame()

# Run PCA
pca <- prcomp(cpm)

# Get percent of variance explained by each PC
PoV <- pca$sdev^2/sum(pca$sdev^2) * 100  
PoV <- round(PoV, digits = 1)

# Merge with sample info, then plot
pca <- pca$rotation |> 
  as.data.frame() |> 
  dplyr::select(PC1, PC2) 
pca$new.id <- row.names(pca)

my_palette <- c( "#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00")
legend.names <- c("Sham_1","Sham_2", "MI_1", "MI_2")


#! Add an arg to specify which samples should have labels
pca <- pca |> left_join(pheno.reorder) |> 
  mutate(gene_treat = factor(paste(genotype, treatment), 
              levels = c("WT Sham", "WT MI", "cmAKO Sham", "cmAKO MI"))) |> 
  filter(type == "whole") 

pca.labels <- pca |> filter(new.id == "WT Sham (4)")
pca.plot <- pca |> 
  ggplot(aes(x = PC1, y = PC2, color = gene_treat)) +
  geom_point(size = 8, color = "black") +
  geom_point(size = 7) +
  geom_label_repel(data = pca.labels, aes(x = PC1, y = PC2),
                   label = pca.labels$new.id, color = "black",
                   alpha = 0.9, fill = "#A6CEE3", 
                   box.padding = 0.5, segment.curvature = -0.2,
                   segment.ncp = 3, segment.angle = 20, force = 10, 
                   max.overlaps = 8, force_pull = 0.001, size = 6,
                   nudge_x = 0.015, nudge_y = 0.02) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_x_continuous(expand = expansion(mult = 0.3), name = paste0("PC1", " (", PoV[1], " % of total variance)")) +
  scale_y_continuous(expand = expansion(mult = 0.3), name = paste0("PC2", " (", PoV[2], " % of total variance)")) +
  theme(axis.text.x = element_text(vjust = 0.5),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.justification = c("center", "center"),
        legend.box.just = "center",
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'),
        plot.margin = unit(c(1,1,1,1), units = "cm"),
        text = element_text(size = 25)) 

pca.plot
# Save plot to results 
png(file = "results/10_plot_de/pca.png",
    width = 10, 
    height = 8,
    units = "in",
    res = 800)

pca.plot

dev.off()

