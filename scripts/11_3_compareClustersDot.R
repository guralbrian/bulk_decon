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


# Highlight points that became significant for either cell type
# For each ID, tag it as sig for either cell type
cell.go <- go |> 
  select(ID, contrast, qscore.adj) |> 
  pivot_wider(names_from = contrast, values_from = qscore.adj) 
cell.go[is.na(cell.go)] <- 0 
cell.go <- cell.go |> 
  mutate(cell.type.sig = case_when(
    clr.Cardiomyocytes > 1.3  & clr.Fibroblast < 1.3 ~ "Cardiomyocytes",
    clr.Cardiomyocytes < 1.3  & clr.Fibroblast > 1.3 ~ "Fibroblasts",
    clr.Cardiomyocytes > 1.3  & clr.Fibroblast > 1.3 ~ "Both",
    clr.Cardiomyocytes < 1.3  & clr.Fibroblast < 1.3 ~ "Neither"
  )) |> 
  select(ID, cell.type.sig) |> 
  right_join(go)

p.dot <- cell.go |> 
  filter(!(contrast %in% c("clr.Fibroblast", "clr.Cardiomyocytes", "genotype_cmAKO_vs_WT"))) |> 
  ggplot(aes(x = qscore.unadj, y = qscore.adj, color = cell.type.sig)) +
  geom_jitter(alpha = 0.4, size = 4, width = 0.35, height = 0.35) +
  geom_abline(intercept = 0, slope = 1)+
  facet_wrap(~contrast, 
             labeller = labeller(contrast = c("treatment_MI_vs_Sham" = "Myocardial Infarction", 
                                              "genotype_cmAKO_vs_WT" = "cmAKO",
                                              "treatmentMI.genotypecmAKO" = "cmAKO:MI"))) +
  scale_color_manual(values = c("#F96A5F","#66C2A5", "#6683D4", "#828282"), name = "Variables significant\n      for GO term") +
  labs(x = "Q-score without cell types", 
       y = "Q-score with cell types") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_blank(),
    text = element_text(size = 15)
  ) +
  guides(color = guide_legend(nrow = 2))


if(!dir.exists("results/11_clusterProfiler")){
  dir.create("results/11_clusterProfiler")
}

# Save plot to results 
png(file = paste0("results/11_clusterProfiler/go_contrast_dot.png"),
    width = 6, 
    height =4,
    units = "in",
    res = 600)

p.dot

dev.off()

