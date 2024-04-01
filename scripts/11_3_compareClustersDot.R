# Make delta plots for clusterProfiler results

# Load libs
libs <- c("tidyverse", "clusterProfiler","patchwork","data.table") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

go.adj <- readRDS("data/processed/pathway_genesets/goadjusted_any_p.RDS")
go.unadj <- readRDS("data/processed/pathway_genesets/gounadjusted_any_p.RDS")

go.adj <- lapply(names(go.adj),
                       function(x){
                         data <- go.adj[[x]]@result
                         data <- data |> mutate(qscore = -log(p.adjust, base=10))
                         data$contrast <- x
                         data$model <- "adj"
                         return(data)}) |> 
  rbindlist() |> 
  as.data.frame()

go.unadj <- lapply(names(go.unadj),
                 function(x){
                   data <- go.unadj[[x]]@result
                   data <- data |> mutate(qscore = -log(p.adjust, base=10))
                   data$contrast <- x
                   data$model <- "unadj"
                   return(data)}) |> 
  rbindlist() |> 
  as.data.frame()

go <- full_join(go.adj, go.unadj, suffix = c(".adj", ".unadj"), by = c("ID", "contrast", "Description"))
rm(go.adj, go.unadj)
#go <- readRDS("data/processed/pathway_genesets/go_all_simp_005.RDS")
# Check if GO terms have up and down regulation
# Merge matched datasets in go
# Make into dataframes of result with column for contrast and model version

p.dot <- go |> 
  filter(!(contrast %in% c("clr.Fibroblast", "clr.Cardiomyocytes"))) |> 
  #mutate(sig = case_when())
  ggplot(aes(x = qscore.unadj, y = qscore.adj)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)+
  facet_wrap(~contrast) +
  labs(title = "GO Term Significance") +
  theme_minimal()


if(!dir.exists("results/11_clusterProfiler")){
  dir.create("results/11_clusterProfiler")
}

# Save plot to results 
png(file = paste0("results/11_clusterProfiler/go_contrast_dot.png"),
    width = 8, 
    height =4,
    units = "in",
    res = 300)

p.dot

dev.off()