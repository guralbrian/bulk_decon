# Make venn diagram of DESeq2 results overlaps
libs <- c("tidyverse", "DESeq2", "VennDiagram") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load the results and expression matrix

# Results from model without cell type proportions
unadj.res <- readRDS("data/processed/models/unadjusted_de_interaction.RDS")  
res.types <- resultsNames(unadj.res)[-1]

pullResults <- function(comparison){
  data <- results(unadj.res, name=comparison) |> 
  as.data.frame() 
  data <- data[,c("log2FoldChange", "padj")]
  data$type <- comparison
  data$gene <- row.names(data)
  row.names(data) <- NULL
  data}

res.list <- lapply(res.types, pullResults)
names(res.list) <- res.types

# merge the list into a dataframe
res.df <- bind_rows(res.list)

deg.list <- res.df %>%
  filter(padj < 0.05) %>%
  split(.$type) %>%
  lapply(function(df) df$gene)

# Generate the Venn diagram for the contrasts (up to 5)
group.names <- c("Genotype", "MI Treatment", "Interaction of Genotype\nand Treatment")
venn.plot <- venn.diagram(
  x = deg.list[1:length(deg.list)], 
  category.names = group.names,
  filename = NULL, 
  output = TRUE
)

# Plot the Venn diagram
grid.draw(venn.plot)

