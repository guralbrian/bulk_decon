# Make delta plots for clusterProfiler results

# Load libs
libs <- c("tidyverse", "clusterProfiler","patchwork","data.table") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Check if the simplified GO output exists, make it if not
if(!file.exists("data/processed/pathway_genesets/go_all_simp_005.RDS")){
# Load data
go.adj <- readRDS("data/processed/pathway_genesets/goadjusted_005.RDS")
go.unadj <- readRDS("data/processed/pathway_genesets/gounadjusted_005.RDS")

# Make names unique and join into one obj
names(go.adj) <- paste(names(go.adj), "adj", sep = "_")
names(go.unadj) <- paste(names(go.unadj), "unadj", sep = "_")
go <- c(go.adj, go.unadj)
rm(go.adj, go.unadj)

# Simplify terms by semantic similarity
go.simp <- lapply(go, clusterProfiler::simplify, cutoff = 0.6)
saveRDS(go.simp, "data/processed/pathway_genesets/go_all_simp_005.RDS")
}else{
go.simp <- readRDS(go.simp, "data/processed/pathway_genesets/go_all_simp_005.RDS")
}

# Check if GO terms have up and down regulation
# Merge matched datasets in go
# Make into dataframes of result with column for contrast and model version
names.treament <- names(go.simp)[names(go.simp) |> str_detect(pattern = "treatmentMI.genotypecmAKO")]
## tHis could be a function applied to all names
go.treatment <- lapply(names.treament,
                       function(x){
                       data <- go.simp[[x]]@result
                       data$contrast <- x
                       return(data)}) |> 
                rbindlist() |> 
                as.data.frame()

common.go.terms <- go.treatment[duplicated(go.treatment$ID),"ID"]
go.treatment <- go.treatment |> filter(ID %in% common.go.terms)


plotGO <- function(x){
  top.change <- x |>  
    mutate(qscore = -log(p.adjust, base=10)) |> 
    group_by(ID) |> 
    mutate(min = min(qscore),
           max = max(qscore),
           range = max - min) |> 
    ungroup() |> 
    filter(contrast == unique(x$contrast)[[1]]) |> 
    arrange(desc(range)) |> 
    slice_head(n = 15) |> 
    pull(ID)
  top.terms <- x |>  
    filter(contrast == unique(x$contrast)[[1]]) |> 
    arrange(p.adjust) |> 
    slice_head(n = 30) |> 
    pull(ID)
  df <- x |> 
    mutate(qscore = -log(p.adjust, base=10),
           desc.wrap = str_to_title(Description) |> 
             str_wrap(width = 18) |> 
             factor()) |> 
    filter(ID %in% top.change) |>
    arrange(desc(qscore))
   
  df$desc.wrap <- factor(df$desc.wrap, levels = rev(unique(df$desc.wrap)))
  
  p.ego <- ggplot(df, aes(x = desc.wrap, y = qscore, fill = contrast, group = contrast)) +
    geom_bar(position = "dodge", stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("#d16a0a", "#30acd1"), labels = c("Post-Adjustment", "Pre-Adjustment")) +
    theme(
      #axis.text.y = element_text(hjust = 0.5, size = 12),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "darkgray"),
      legend.position = "bottom",
      #legend.key.size = unit(1, "in"),
      #title = element_text(size = 18)
      text = element_text(size = 12)
    ) +
    labs(fill = "Adjusted\np-value",
         y = "Q Score")

}

p.treat <- plotGO(go.treatment)

if(!dir.exists("results/11_clusterProfiler")){
  dir.create("results/11_clusterProfiler")
}

# Save plot to results 
png(file = paste0("results/11_clusterProfiler/go_inter_contrast.png"),
    width = 8, 
    height =12,
    units = "in",
    res = 100)

p.treat

dev.off()

# Sankey Diagram

# Library
library(networkD3)
library(dplyr)

# Merge go datasets
# Make into dataframes of result with column for contrast and model version

go.merged <- lapply(names(go.simp),
                       function(x){
                         data <- go.simp[[x]]@result
                         data$contrast <- x
                         return(data)}) |> 
                    rbindlist() |> 
                    as.data.frame() |> 
                    select(ID, Description, p.adjust, contrast)

# Make a before and after column from adj/unadj
go.merged <- go.merged |> 
  mutate(
    model = case_when(
      str_detect(contrast, "_unadj") ~ "unadjust",
      .default = "adjusted"
    ),
    contrast = case_when(
      str_detect(contrast, "treatment_MI_vs_Sham") ~ "MI",
      str_detect(contrast, "genotype_cmAKO_vs_WT") ~ "cmAKO",
      str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "MI:cmAKO",
      str_detect(contrast, "Fibroblast") ~ "Fibroblast",
      str_detect(contrast, "Cardio") ~ "Cardiomyocytes"
    )
  )


combinations <- expand.grid(unique(go.merged$model), unique(go.merged$contrast))
combinations <- combinations[c(1:5,7,9,10),]
combination.list <-c()
# Make a loop to pull the IDs from each combination
for(i in 1:nrow(combinations)){
  combination.list[[i]] <- go.merged |> filter(model ==  combinations[i,1] & contrast == combinations[i,2]) |> pull(ID) 
}


names(combination.list) <- paste(combinations$Var1, combinations$Var2, sep = "_")

comb.table <- crossprod(table(stack(combination.list))) |> as.data.frame()
comb.table$source <- row.names(comb.table)

comb.table <-pivot_longer(comb.table, cols = names(combination.list), values_to = "value", names_to = "target")

comb.table <- comb.table |> mutate(
  filt = case_when(
      str_detect(source, "unadjust") & str_detect(target, "adjusted") ~ TRUE,
      .default = FALSE)) 
comb.table <- comb.table |> filter(filt == TRUE)
#comb.table$source <- lapply(str_split(comb.table$source, "_"), "[[", 2) |> unlist()
#comb.table$target <- lapply(str_split(comb.table$target, "_"), "[[", 2) |> unlist()

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(comb.table$source), 
         as.character(comb.table$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
comb.table$IDsource <- match(comb.table$source, nodes$name)-1 
comb.table$IDtarget <- match(comb.table$target, nodes$name)-1

# Make the Network
p <- sankeyNetwork(Links = comb.table, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)
p
