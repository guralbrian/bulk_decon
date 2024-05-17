# Make tables for GO terms 
# Load libs
libs <- c("tidyverse", "clusterProfiler", "patchwork", "data.table") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

n.terms <- 5
go.adj <- readRDS("data/processed/pathway_genesets/goadjusted_005.RDS")
go.unadj <- readRDS("data/processed/pathway_genesets/gounadjusted_005.RDS")


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

# Make table for each comparison
# Columns needed: GO ID, description, P-value & Qscore, before /after
options(digits=3)

## Make a table for each:
# GO terms before adj, most significant, no unsig
# Top terms after adj for each term (5 terms)
# Most increased and decreased terms

#### Most Decreased qscores ####
# make plot for factor comparisons

# Subset to variables of interest, get most decreased in sig
decreas.dat <- go |>
  filter(qscore.adj > 1.3 | qscore.unadj > 1.3) |> 
  select(Description, p.adjust.unadj, p.adjust.adj, qscore.unadj, qscore.adj,  contrast) |> 
  filter(contrast %in% 
           c("treatment_MI_vs_Sham")) |> 
             #"genotype_cmAKO_vs_WT", 
            # c("treatmentMI.genotypecmAKO")) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description) |> str_wrap(width = 25),
         delta_q = qscore.adj- qscore.unadj,
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI"
         )) |> 
  arrange(desc(qscore.adj)) |> 
  slice_head(n = 10) 

# get top terms in order of decreased sig
decrease.order <- decreas.dat |>
  arrange(qscore.adj) |> 
  pull(Description)

# Subset to Descriptions to match top plot and get cell type values
decrease.cell <- go |>
  #filter(qscore.adj > 1.3 | qscore.unadj > 1.3) |> 
  select(Description, p.adjust.adj , qscore.adj,  contrast) |> 
  filter(contrast %in% 
           c("clr.Fibroblast", 
             "clr.Cardiomyocytes")) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description) |> str_wrap(width = 25),
         contrast = case_when(
           str_detect(contrast, "clr.Cardiomyocytes") ~ "Cardiomyocytes",
           str_detect(contrast, "clr.Fibroblast") ~ "Fibroblasts"
         )) |> 
  filter(Description %in% decrease.order) 

# get top of axis range for each plot, take top

q.max.cell <- decrease.cell |> 
  ungroup() |> 
  arrange(desc(qscore.adj)) |> 
  slice_head(n=1) |> 
  pull(qscore.adj)

q.max.var <- decreas.dat |> 
  pivot_longer(cols = c(qscore.unadj, qscore.adj)) |> 
  ungroup() |> 
  arrange(desc(value)) |> 
  slice_head(n=1) |> 
  pull(value)

xlimit <- max(q.max.cell, q.max.var)
# Barplot
p.de.treat <- decreas.dat |> 
  pivot_longer(cols = c(qscore.unadj, qscore.adj)) |> 
  mutate(name = case_when(
    str_detect(name, "qscore.adj") ~ "Post-adj",
    str_detect(name, "qscore.unadj") ~ "Pre-adj",
  )) |> 
  ggplot(aes(x = factor(Description, levels = as.character(decrease.order)), y = value, fill = name)) + 
  geom_bar(stat = "summary", position = position_dodge(0.9),  fun = mean,width = 0.9,  color = "black", alpha = 0.8) +
  #facet_wrap(~contrast, scales = "free", ncol =1) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    legend.position = "bottom",
    axis.text.y.left = element_text(margin = margin(0, 1, 0, 1), hjust = 0.5)
  )+
  scale_fill_manual(values = c("#A6D854", "#D95F02")) +
  scale_y_continuous(limits = c(0, xlimit)) +
  labs(y = "Q-value (-log10[p-value])", fill = "")



p.de.cell <- decrease.cell |> 
  ggplot(aes(x = factor(Description, levels = as.character(decrease.order)), 
             y = qscore.adj, 
             fill = factor(contrast, levels = c("Fibroblasts", "Cardiomyocytes")))) + 
  geom_bar(stat = "summary", position = position_dodge(0.9),  fun = mean,width = 0.9,  color = "black", alpha = 0.8) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_blank(), 
    axis.ticks = element_blank(),
    plot.margin = unit(c(1,0,1,-1), "mm") 
  ) +
  scale_fill_manual(values = c("#8DA0CB","#66C2A5")) +
  scale_y_continuous(limits = c(0, xlimit)) +
  scale_y_reverse() +
  labs(fill = "")



png(file = "results/misc/goBarInfarctTopSig.png",
    width = 9, 
    height = 6,
    units = "in",
    res = 600)

wrap_plots(p.de.cell, p.de.treat)


dev.off()


# Need to:
# Change colors for cells vs comparison plot
# Make qscore axis consistant
# Center middle text 
# repeat for each comparison 
# repeat for:
## Most increased


##### Plots for pre-adj
# Subset to variables of interest, get most decreased in sig
pre.adj <- go |>
  select(Description, qscore.unadj,  contrast) |> 
  filter(contrast %in% 
           c("treatment_MI_vs_Sham", "treatmentMI.genotypecmAKO")) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description) |> str_wrap(width = 25),
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI"
         )) |> 
  arrange(desc(qscore.unadj)) |> 
  slice_head(n = 6) 

# get top terms in order of decreased sig
decrease.order <- pre.adj |>
  arrange(qscore.unadj) |> 
  pull(Description)


# Barplot
p.preadj <- pre.adj |>
  ggplot(aes(x = factor(Description, levels = as.character(decrease.order)), y = qscore.unadj, fill = contrast)) + 
  geom_bar(stat = "summary", position = position_dodge(0.9),  fun = mean,width = 0.9,  color = "black", alpha = 0.8) +
  facet_wrap(~contrast, scales = "free", ncol =2) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.y.left = element_text(margin = margin(0, 1, 0, 1), hjust = 0.5),
    strip.text = element_text(size = 12)
  )+
  scale_fill_manual(values = c("#A6D854", "#D95F02")) +
  labs(y = "Q-value (-log10[p-value])", fill = "")


png(file = "results/misc/goBarPreAdjTopSig.png",
    width = 8, 
    height = 4,
    units = "in",
    res = 600)
p.preadj

dev.off()
