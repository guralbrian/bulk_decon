# Upset plot of DESeq2 results
libs <- c("tidyverse", "DESeq2","ComplexUpset", "patchwork") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Results from model without cell type proportions
adj.res <- readRDS("data/processed/models/adjusted_de_interaction.RDS")  
res.types <- resultsNames(adj.res)[-1]

# Prepare a data frame to collect significant gene information
significant_genes_df <- data.frame(gene = character(), condition = character(), regulation = character(), stringsAsFactors = FALSE)

# Loop through each type of result to get significant genes, split by up or down regulation
for (type in res.types) {
  res <- results(adj.res, name = type)
  
  # Filter for significant genes
  sig <- subset(res, padj < 0.05)
  
  # Determine up or down regulation
  sig$regulation <- ifelse(sig$log2FoldChange > 0, "Up", "Down")
  
  # Add to the data frame
  tmp_df <- data.frame(gene = rownames(sig), condition = type, regulation = sig$regulation, stringsAsFactors = FALSE)
  significant_genes_df <- rbind(significant_genes_df, tmp_df)
}

# Create a binary matrix for the UpSet plot
# Each row represents a gene, and each column represents a condition (combination of type and regulation)
down.upset <- significant_genes_df %>%
  filter(regulation == "Down")|> 
  select(-regulation) |> 
  table() |> 
  as.data.frame.matrix()

up.upset <- significant_genes_df %>%
  filter(regulation == "Up")|> 
  select(-regulation) |> 
  table() |> 
  as.data.frame.matrix()

colnames(up.upset) <- colnames(down.upset) <- c("CMs", "FBs", "cmAKO", "MI", "cmAKO X MI")


presence = ComplexUpset:::get_mode_presence('exclusive_intersection')

summarise_values = function(df) {
  aggregate(
    as.formula(paste0(presence, '~ intersection')),
    df,
    FUN=sum
  )
}

p.upset.down <- down.upset |> 
  ComplexUpset::upset(intersect = colnames(down.upset),
                      base_annotations=list(
                        'log10(intersection size)'=(
                          ggplot()
                          + geom_bar(data=summarise_values, stat='identity', aes(y=!!presence), 
                                     color = "black", fill = "#d15469") 
                          + scale_y_continuous(trans='log10')
                          + labs(title = "Downregulated")
                          + theme(text = element_text(size = 10),
                                  plot.margin = unit(c(0,0,0,0), unit = "cm"),
                                  axis.title = element_blank(),
                                  axis.text.y = element_blank())
                          
                        )
                      ),
                      labeller=ggplot2::as_labeller(c(
                        'CMs' = "",
                        'FBs' = "",
                        'cmAKO'='',
                        'MI'='',
                        'cmAKO X MI' = ""
                      )),
                      min_size=5,
                      width_ratio=0.1,
                      set_sizes = F
  )

p.upset.up <- up.upset |> 
  ComplexUpset::upset(intersect = colnames(up.upset),
                      base_annotations=list(
                        'log10(intersection size)'=(
                          ggplot()
                          + geom_bar(data=summarise_values, stat='identity', aes(y=!!presence), 
                                     color = "black", fill = "lightblue") 
                          + scale_y_continuous(trans='log10') 
                          + labs(title = "Upregulated", y = 'Sig. Genes')
                          + theme(text = element_text(size = 10),
                                  plot.margin = unit(c(0,0,0,0), unit = "cm"),
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_text(vjust = -20))
                        )
                      ),
                      wrap = F,
                      min_size=5,
                      width_ratio=0.1,
                      set_sizes = F
  )


# Save plot to results 

if(!dir.exists("results/10_plot_de")){
  dir.create("results/10_plot_de")
}

png(file = "results/10_plot_de/upset_adj.png",
    width = 6, 
    height = 3,
    units = "in",
    res = 600)

wrap_plots(p.upset.up, p.upset.down)

dev.off()

