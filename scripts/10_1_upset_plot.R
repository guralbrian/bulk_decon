# Upset plot of DESeq2 results
libs <- c("tidyverse", "DESeq2","ComplexUpset") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Results from model without cell type proportions
unadj.res <- readRDS("data/processed/models/unadjusted_de_interaction.RDS")  
res.types <- resultsNames(unadj.res)[-1]

# Prepare a data frame to collect significant gene information
significant_genes_df <- data.frame(gene = character(), condition = character(), regulation = character(), stringsAsFactors = FALSE)

# Loop through each type of result to get significant genes, split by up or down regulation
for (type in res.types) {
  res <- results(unadj.res, name = type)
  
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
d.upset <- significant_genes_df %>%
  filter(regulation == "Down")|> 
  select(-regulation) |> 
  table() |> 
  as.data.frame.matrix()


presence = ComplexUpset:::get_mode_presence('exclusive_intersection')

summarise_values = function(df) {
  aggregate(
    as.formula(paste0(presence, '~ intersection')),
    df,
    FUN=sum
  )
}

p.upset <-d.upset |> 
  ComplexUpset::upset(intersect = colnames(d.upset),
                      base_annotations=list(
    'log10(intersection size)'=(
      ggplot()
      + geom_bar(data=summarise_values, stat='identity', aes(y=!!presence))
      + ylab('Significant Genes')
      + scale_y_continuous(trans='log10')
    )
  ),
  min_size=5,
  width_ratio=0.1
  )
p.upset
