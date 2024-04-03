# Make tables for GO terms 
# Load libs
libs <- c("tidyverse", "clusterProfiler", "gt", "data.table") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#go.adj <- readRDS("data/processed/pathway_genesets/goadjusted_any_p.RDS")
#go.unadj <- readRDS("data/processed/pathway_genesets/gounadjusted_any_p.RDS")

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
tab <- go |>
  filter(qscore.adj > 1.3 | qscore.unadj > 1.3) |> 
  filter(contrast %in% 
           c("treatment_MI_vs_Sham", 
             "genotype_cmAKO_vs_WT", 
             "treatmentMI.genotypecmAKO")) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description),
         delta_q = qscore.adj- qscore.unadj,
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI"
         )) |> 
  arrange(delta_q) |> 
  slice_head(n = 5) |> 
  select(Description, p.adjust.unadj, p.adjust.adj, qscore.unadj, qscore.adj, delta_q, contrast) |> 
  gt(rowname_col = "Description",
     groupname_col = "contrast") |>
  tab_spanner(
    label = "Cell-type adjusted",
    columns = c(p.adjust.adj, qscore.adj)
  ) |> 
  tab_spanner(
    label = "Cell-type not considered",
    columns = c(p.adjust.unadj, qscore.unadj)
  ) |> 
  cols_merge(
    columns = c(p.adjust.unadj, qscore.unadj),
    pattern = "{1} ({2})"
  ) |>
  cols_merge(
    columns = c(p.adjust.adj, qscore.adj),
    pattern = "{1} ({2})"
  ) |> 
  cols_label(
    p.adjust.unadj = "p-value (Q-score)",
    p.adjust.adj = "p-value (Q-score)",
    delta_q = paste0(html("\u394"), "Q-score")
  ) |> 
  cols_align(
    align = "center",
    columns = c(p.adjust.unadj, p.adjust.unadj, p.adjust.adj, p.adjust.adj)
  ) |> 
  tab_header(
    title = md("GO terms associated cmAKO, MI, and their interaction"),
    subtitle = "By significance lost after cell-type-adjustment"
  ) |>  
  opt_row_striping() |> 
  cols_width(Description ~ 200) |>  
  tab_style(
    style = list(
      align = "center",
      cell_fill("grey"),
      cell_text(color = "black", weight = "bold")),
    locations = cells_row_groups())

file = "results/11_clusterProfiler/tables/lost_q_table"
gtsave(tab, paste0(file, ".html"))

webshot::webshot(url = paste0(file, ".html"), file = paste0(file, ".png"), 
                 vwidth = 700, vheight = 1200)

#### Most Decreased qscores ####
tab <- go |>
  filter(qscore.adj > 1.3 | qscore.unadj > 1.3) |> 
  filter(contrast %in% 
           c("treatment_MI_vs_Sham", 
             "genotype_cmAKO_vs_WT", 
             "treatmentMI.genotypecmAKO")) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description),
         delta_q = qscore.adj- qscore.unadj,
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI"
         )) |> 
  arrange(desc(delta_q)) |> 
  slice_head(n = 5) |> 
  select(Description, p.adjust.unadj, p.adjust.adj, qscore.unadj, qscore.adj, delta_q, contrast) |> 
  gt(rowname_col = "Description",
     groupname_col = "contrast") |>
  tab_spanner(
    label = "Cell-type adjusted",
    columns = c(p.adjust.adj, qscore.adj)
  ) |> 
  tab_spanner(
    label = "Cell-type not considered",
    columns = c(p.adjust.unadj, qscore.unadj)
  ) |> 
  cols_merge(
    columns = c(p.adjust.unadj, qscore.unadj),
    pattern = "{1} ({2})"
  ) |>
  cols_merge(
    columns = c(p.adjust.adj, qscore.adj),
    pattern = "{1} ({2})"
  ) |> 
  cols_label(
    p.adjust.unadj = "p-value (Q-score)",
    p.adjust.adj = "p-value (Q-score)",
    delta_q = paste0(html("\u394"), "Q-score")
  ) |> 
  cols_align(
    align = "center",
    columns = c(p.adjust.unadj, p.adjust.unadj, p.adjust.adj, p.adjust.adj)
  ) |> 
  tab_header(
    title = md("GO terms associated cmAKO, MI, and their interaction"),
    subtitle = "By significance gained after cell-type-adjustment"
  ) |>  
  opt_row_striping() |> 
  cols_width(Description ~ 200) |>  
  tab_style(
    style = list(
      align = "center",
      cell_fill("grey"),
      cell_text(color = "black", weight = "bold")),
    locations = cells_row_groups()) 

file = "results/11_clusterProfiler/tables/upped_q_table"
gtsave(tab, paste0(file, ".html"))

webshot::webshot(url = paste0(file, ".html"), file = paste0(file, ".png"), 
                 vwidth = 700, vheight = 1200)
#### Top terms, unadj ####
tab <- go |> 
  filter(contrast %in% 
           c("treatment_MI_vs_Sham", 
             "genotype_cmAKO_vs_WT", 
             "treatmentMI.genotypecmAKO")) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description),
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI"
         )) |> 
  filter(qscore.unadj >= 1.3) |> 
  arrange(desc(qscore.unadj)) |> 
  slice_head(n = 5) |> 
  select(Description, p.adjust.unadj, qscore.unadj,  contrast) |> 
  gt(rowname_col = "Description",
     groupname_col = "contrast") |>
  cols_merge(
    columns = c(p.adjust.unadj, qscore.unadj),
    pattern = "{1} ({2})"
  ) |>
  cols_label(
    p.adjust.unadj = "p-value (Q-score)"
  ) |> 
  cols_align(
    align = "center",
    columns = c(p.adjust.unadj, qscore.unadj)
  ) |> 
  tab_header(
    title = md("GO terms associated cmAKO, MI, and their interaction"),
    subtitle = "Top five terms by significance in overrepresentation testing"
  ) |>  
  opt_row_striping() |> 
  cols_width(Description ~ 300) |> 
  tab_style(
    style = list(
      align = "center",
      cell_fill("grey"),
      cell_text(color = "black", weight = "bold")),
    locations = cells_row_groups())

file = "results/11_clusterProfiler/tables/unadj_top_table"
gtsave(tab, paste0(file, ".html"))
webshot::webshot(url = paste0(file, ".html"), file = paste0(file, ".png"), 
                 vwidth = 600, vheight = 800)

#### Top terms, adj ####
tab <- go |> 
  filter(qscore.adj >= 1.3) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description),
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI",
           str_detect(contrast, "clr.Cardiomyocytes") ~ "Cardiomyocytes",
           str_detect(contrast, "clr.Fibroblast") ~ "Fibroblasts"
         )) |>
  arrange(desc(qscore.unadj)) |> 
  slice_head(n = 5) |> 
  select(Description, p.adjust.adj, qscore.adj,  contrast) |> 
  gt(rowname_col = "Description",
     groupname_col = "contrast") |>
  cols_merge(
    columns = c(p.adjust.adj, qscore.adj),
    pattern = "{1} ({2})"
  ) |>
  cols_label(
    p.adjust.adj = "p-value (Q-score)"
  ) |> 
  cols_align(
    align = "center",
    columns = c(p.adjust.adj, qscore.adj)
  ) |> 
  tab_header(
    title = md("GO terms associated cmAKO, MI, and their interaction"),
    subtitle = "Top five terms by significance after adjusting for cell-types"
  ) |>  
  opt_row_striping() |> 
  cols_width(Description ~ 300) |> 
  tab_style(
    style = list(
      align = "center",
      cell_fill("grey"),
      cell_text(color = "black", weight = "bold")),
    locations = cells_row_groups())

file = "results/11_clusterProfiler/tables/adj_top_table"
gtsave(tab, paste0(file, ".html"))
webshot::webshot(url = paste0(file, ".html"), file = paste0(file, ".png"), 
                 vwidth = 600, vheight = 1200)
