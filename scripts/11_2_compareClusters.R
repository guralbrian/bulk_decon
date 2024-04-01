# Make barplots of top changed clusters in either direction 

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


## Processing for increased significance GO terms

go.plot <- go |> 
  filter(qscore.unadj > 2 | qscore.adj > 2 ) |> 
  filter(!(contrast %in% c("clr.Fibroblast", "clr.Cardiomyocytes"))) |>
  mutate(delta.q = qscore.unadj - qscore.adj,
         desc.wrap = str_to_title(Description) |> 
           str_wrap(width = 18) |> 
           factor()) |> 
  group_by(contrast) |> 
  arrange(delta.q) |> 
  slice_head(n = 6) |> 
  pivot_longer(cols = c("qscore.adj", "qscore.unadj"))
desc.wrap <- go.plot |> 
  filter(name == "qscore.adj") |> 
  arrange(desc(value)) |> 
  pull(desc.wrap)

go.plot$desc.wrap <- factor(go.plot$desc.wrap, levels = rev(unique(desc.wrap)))
p.go.up <- go.plot |> 
  ggplot(aes(x = desc.wrap, y = value, fill = name, group = name)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~contrast, scales = "free") +
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
       y = "Q Score",
       title = "Top GO terms by gained significance")


wrap_plots(p.treat.unadj, p.treat.adj)
if(!dir.exists("results/11_clusterProfiler")){
  dir.create("results/11_clusterProfiler")
}

# Save plot to results 
png(file = paste0("results/11_clusterProfiler/go_increased_sig_bar.png"),
    width = 10, 
    height =6,
    units = "in",
    res = 300)

p.go.up

dev.off()


# Scripts for decreased 

go.plot <- go |> 
  filter(qscore.unadj > 2 | qscore.adj > 2 ) |> 
  filter(!(contrast %in% c("clr.Fibroblast", "clr.Cardiomyocytes"))) |>
  mutate(delta.q = qscore.unadj - qscore.adj,
         desc.wrap = str_to_title(Description) |> 
           str_wrap(width = 18) |> 
           factor()) |> 
  group_by(contrast) |> 
  arrange(desc(delta.q)) |> 
  slice_head(n = 6) |> 
  pivot_longer(cols = c("qscore.adj", "qscore.unadj"))
desc.wrap <- go.plot |> 
  filter(name == "qscore.unadj") |> 
  arrange(desc(value)) |> 
  pull(desc.wrap)

go.plot$desc.wrap <- factor(go.plot$desc.wrap, levels = rev(unique(desc.wrap)))
p.go.down <- go.plot |> 
  ggplot(aes(x = desc.wrap, y = value, fill = name, group = name)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~contrast, scales = "free") +
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
       y = "Q Score",
       title = "Top GO terms by lost significance")


if(!dir.exists("results/11_clusterProfiler")){
  dir.create("results/11_clusterProfiler")
}

# Save plot to results 
png(file = paste0("results/11_clusterProfiler/go_decreased_sig_bar.png"),
    width = 10, 
    height =6,
    units = "in",
    res = 300)

p.go.down

dev.off()
