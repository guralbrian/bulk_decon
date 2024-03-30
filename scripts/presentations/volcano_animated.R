# Visualize DE 
libs <- c("tidyverse","ggrepel", "DESeq2","patchwork", "gganimate", "data.table") # list libraries here
lapply(libs, library, character.only = T, quietly = T)
rm(libs)

# Load the results 
# Load the results and expression matrices
res.unadj <- readRDS(paste0("data/processed/models/unadjusted_de_interaction.RDS")) 
res.adj <- readRDS(paste0("data/processed/models/adjusted_de_interaction.RDS")) 

# Function to pull results from each RDS obj
processDE <- function(des.res){
  names <- resultsNames(des.res)[-1]
  # make a list of each result
  des.res <- lapply(names, function(x){
    results(des.res, name=x) |> 
      as.data.frame() |> 
      rownames_to_column(var = "gene")
  })
  names(des.res) <- names
  return(des.res)
}

# Get list of DE results
res.unadj <- processDE(res.unadj)
res.adj <- processDE(res.adj)

# Make a dataframe of DEGs that includes pre- and post-adjustment values

overlaps <- names(res.adj)[names(res.adj) %in% names(res.unadj)]

# Add column to indicate name
adj.named <- lapply(res.adj, function(x){
      x <- as.data.frame(x)
      x$comparison <- "adjusted"
      return(x)})
unadj.named <- lapply(res.unadj, function(x){
  x <- as.data.frame(x)
  x$comparison <- "unadjusted"
  return(x)})

# Merge the two with a loop (oh loops..)
bound.res <- c()
for(i in overlaps){
  print(i)
  bound.res[[i]] <- rbind(unadj.named[[i]], adj.named[[i]])
  bound.res[[i]][["variable"]] <- i
}

# Collapse into a single dataframe
res <- rbindlist(bound.res)
write.csv(res, "data/processed/models/all_deseq.csv")
genes <- c("Errfi1", "Egr1", "Pik3r1", "Zbtb16")

top <- res |> filter(gene %in% genes & variable == "treatmentMI.genotypecmAKO")
# Find top DE genes
res <- res |> 
  mutate(target = case_when(
    gene %in% genes ~ "1",
    .default = "0"),
         variable = as.factor(variable),
         comparison = as.factor(comparison)) |> 
  filter(variable == "treatmentMI.genotypecmAKO")


plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = target)) +
  geom_point(alpha = 0.8, size = 8) + 
  scale_color_manual(values = c("#999999", "#ed9209")) +
  geom_text(data = top, aes(label = gene), size = 14, color = "black", nudge_y = 0.7) + # label most sig genes 
  #facet_wrap(vars(variable)) +
  theme_minimal() +
  xlim(-5, 5) +
  labs(x = "log2(Fold Change)", y = "-log10(adjusted p-value)", 
       color = "p < 0.05 and\nFold Change > 1.5") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "solid") + # add a line for p=0.05
  theme(legend.position = "none",
        axis.text = element_text(color = "black", size = 30),
        axis.ticks = element_blank(),
        title = element_text(size = 35),
        axis.title = element_text(color = "black", size = 30),
        panel.background = element_rect(color="black"),
        plot.background = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot
options(gganimate.dev_args = list(width = 1600, height = 900, fps = 30))
p.ani <- plot + transition_states(comparison, transition_length = 5) +
  labs(title = "Cell types {closest_state}")
anim_save("results/mohlke/animation_test.gif", p.ani, width = 1000, height = 900, nframes = 300, fps = 50)
