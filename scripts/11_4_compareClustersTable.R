# Make tables for GO terms 
# Load libs
libs <- c("tidyverse", "clusterProfiler", "gt", "data.table") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

n.terms <- 5
go.adj <- readRDS("data/processed/pathway_genesets/goadjusted_005.RDS")
go.unadj <- readRDS("data/processed/pathway_genesets/gounadjusted_005.RDS")

# Simplify terms


go.adj <- lapply(names(go.adj),
                 function(x){
                   data <- go.adj[[x]]@result
                   data <- data |> mutate(qscore = -log(p.adjust, base=10))
                   data$contrast <- x
                   return(data)}) |> 
  rbindlist() |> 
  as.data.frame()

go.unadj <- lapply(names(go.unadj),
                   function(x){
                     data <- go.unadj[[x]]@result
                     data <- data |> mutate(qscore = -log(p.adjust, base=10))
                     data$contrast <- x
                     return(data)}) |> 
  rbindlist() |> 
  as.data.frame()

go <- full_join(go.adj, go.unadj, suffix = c(".adj", ".unadj"), by = c("ID", "contrast", "Description"))
#rm(go.adj, go.unadj)


# Make a supplement dataset
write.csv(go, "results/supp_data/go_terms.csv")
options(digits=3)

## Make a table for each:
# GO terms before adj, most significant, no unsig
# Top terms after adj for each term (5 terms)
# Most increased and decreased terms


#### Top terms, unadj ####
tab <- go |> 
  filter(contrast %in% 
           c("treatment_MI_vs_Sham", 
             "treatmentMI.genotypecmAKO")) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description),
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO X Myocardial Infarction"
         )) |> 
  filter(qscore.unadj >= 1.3) |> 
  group_by(contrast, qscore.unadj) |> # This removes terms with identical p-values
  slice_head(n = 1) |>  # They're effectively just duplicated terms
  group_by(contrast) |> 
  arrange(desc(qscore.unadj)) |> 
  slice_head(n = n.terms) |> 
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
  opt_row_striping() |> 
  cols_width(Description ~ 260) |> 
  tab_style(
    style = list(
      align = "center",
      cell_fill("grey"),
      cell_text(color = "black", weight = "bold")),
    locations = cells_row_groups())

file = "results/11_clusterProfiler/tables/unadj_top_table"
gtsave(tab, paste0(file, ".html"))
webshot::webshot(url = paste0(file, ".html"), file = paste0(file, ".png"), 
                 vwidth = 450, vheight = n.terms*110,  zoom = 3)

#### Top terms, adj ####
tab <- go |> 
  filter(qscore.adj >= 1.3 & contrast != "genotype_cmAKO_vs_WT") |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description),
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI",
           str_detect(contrast, "clr.Cardiomyocytes") ~ "Cardiomyocytes",
           str_detect(contrast, "clr.Fibroblast") ~ "Fibroblasts"
         )) |>
  filter(qscore.adj >= 1.3) |> 
  group_by(contrast, qscore.adj) |> # This removes terms with identical p-values
  slice_head(n = 1) |>  # They're effectively just duplicated terms
  group_by(contrast) |>
  arrange(desc(qscore.adj)) |> 
  slice_head(n = n.terms) |> 
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
                 vwidth = 600, vheight = n.terms*210, zoom = 3)

#### Top terms, adj MI and MI:cmAKO ####
tab <- go |> 
  filter(qscore.adj >= 1.3 & contrast %in% c("treatment_MI_vs_Sham", 
                                             "treatmentMI.genotypecmAKO")) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description),
         delta_q = qscore.adj- qscore.unadj,
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO X Myocardial Infarction"
         )) |>
  group_by(contrast, qscore.adj) |> # This removes terms with identical p-values
  slice_head(n = 1) |>  # They're effectively just duplicated terms
  group_by(contrast) |> 
  arrange(desc(qscore.adj)) |> 
  slice_head(n = n.terms) |> 
  select(Description, p.adjust.unadj, p.adjust.adj, qscore.unadj, qscore.adj, delta_q, contrast) |> 
  gt(rowname_col = "Description",
     groupname_col = "contrast") |>
  tab_spanner(
    label = "Cell-type adjusted",
    columns = c(p.adjust.adj, qscore.adj)
  ) |> 
  tab_spanner(
    label = "Unadjusted",
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
  opt_row_striping() |> 
  cols_width(Description ~ 230) |> 
  tab_style(
    style = list(
      align = "center",
      cell_fill("grey"),
      cell_text(color = "black", weight = "bold")),
    locations = cells_row_groups())

file = "results/11_clusterProfiler/tables/adj_top_treat_gene_table"
gtsave(tab, paste0(file, ".html"))
webshot::webshot(url = paste0(file, ".html"), file = paste0(file, ".png"), 
                 vwidth = 650, vheight = n.terms*110, zoom = 3)

#### Top terms, adj cell types ####
n.terms <- 5
# Get the top 5 terms for each cell type
top.terms <- go |> 
  filter(contrast %in% c("clr.Cardiomyocytes", 
                         "clr.Fibroblast")) |> 
  group_by(contrast, qscore.adj) |> # This removes terms with identical p-values
  slice_head(n = 1) |>
  group_by(contrast) |> 
  arrange(desc(qscore.adj)) |> 
  slice_head(n = n.terms) |> 
  pull(ID)

# find the strongest unadj association of each
top.go <- go |> 
  filter(ID %in% top.terms &
         contrast %in% c("treatment_MI_vs_Sham", 
                         "treatmentMI.genotypecmAKO",
                         "cmAKO_vs_WT")) |> 
  group_by(Description) |> 
  arrange(desc(qscore.unadj)) |> 
  slice_head(n = 1) |> 
  mutate(orig.contrast = contrast,
         orig.p = p.adjust.unadj,
         orig.q = qscore.unadj,
         orig.contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI")) |> 
  select(Description, orig.contrast, orig.p, orig.q)
tab <- go |> 
  filter(contrast %in% c("clr.Cardiomyocytes", 
                         "clr.Fibroblast")) |> 
  group_by(contrast, qscore.adj) |> # This removes terms with identical p-values
  slice_head(n = 1) |>  # They're effectively just duplicated terms
  group_by(contrast) |> 
  arrange(desc(qscore.adj)) |> 
  slice_head(n = n.terms) |> 
  left_join(top.go, by = "Description") |> 
  select(Description, contrast, p.adjust.adj, qscore.adj, orig.p, orig.q, orig.contrast) |> 
  mutate(Description = stringr::str_to_title(Description),
         contrast = case_when(
              str_detect(contrast, "clr.Cardiomyocytes") ~ "Cardiomyocytes",
              str_detect(contrast, "clr.Fibroblast") ~ "Fibroblasts")) |> 
  gt(rowname_col = "Description",
     groupname_col = "contrast") |>
  tab_spanner(
    label = "Cell type assc.",
    columns = c(p.adjust.adj, qscore.adj)
  ) |> 
  tab_spanner(
    label = "Pre-adjustment assc.",
    columns = c(orig.p, orig.q, orig.contrast)
  ) |>
  cols_merge(
    columns = c(p.adjust.adj, qscore.adj),
    pattern = "{1} ({2})") |> 
  cols_merge(
      columns = c(orig.p, orig.q),
      pattern = "{1} ({2})"
    ) |> 
  cols_label(
    orig.p = "p-value (Q-score)",
    p.adjust.adj = "p-value (Q-score)",
    orig.contrast = "Assc. Variable"
    ) |> 
  cols_align(
    align = "center",
    columns = c(p.adjust.adj, p.adjust.adj, orig.p, orig.q, orig.contrast)
  )  |>  
  opt_row_striping() |> 
  cols_width(Description ~ 240) |> 
  tab_style(
    style = list(
      align = "center",
      cell_fill("grey"),
      cell_text(color = "black", weight = "bold")),
    locations = cells_row_groups())

file = "results/11_clusterProfiler/tables/adj_top_cell_table"
gtsave(tab, paste0(file, ".html"))
webshot::webshot(url = paste0(file, ".html"), file = paste0(file, ".png"), 
                 vwidth = 740, vheight = n.terms*110, zoom = 3)

#### Supplement tables ####
#### Fig 4 



#### Most Decreased qscores ####
#### Supp 4b ####
n.terms <- 20
tab <- go |>
  filter(qscore.adj > 1.3 | qscore.unadj > 1.3) |> 
  filter(contrast %in% 
           c("treatment_MI_vs_Sham", 
             #"genotype_cmAKO_vs_WT", 
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
  slice_head(n = n.terms) |> 
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
    title = md("GO terms associated with WT and cmAKO responses to MI"),
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

webshot::webshot(url = paste0(file, ".html"), file =  "results/supp_figs/4b_lost_sig_go.png", 
                 vwidth = 700, vheight = n.terms*110, zoom = 3)

#### Most Increased qscores 
#### Supp 4c ####
n.terms <- 20
tab <- go |>
  filter(qscore.adj > 1.3 | qscore.unadj > 1.3) |> 
  filter(contrast %in% 
           c("treatment_MI_vs_Sham", 
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
  slice_head(n = n.terms) |> 
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
    title = md("Most significantly increased terms after adjustment")
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

# Save to supplements
webshot::webshot(url = paste0(file, ".html"), file = "results/supp_figs/4c_gained_sig_go.png", 
                 vwidth = 700, vheight = n.terms*110, zoom = 3)

#### Longer lists of top terms
#### Supp 4d ####
n.terms <- 20
tab <- go |> 
  filter(qscore.adj >= 1.3 & contrast %in% c("treatment_MI_vs_Sham", "treatmentMI.genotypecmAKO")) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description),
         delta_q = qscore.adj- qscore.unadj,
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI",
           str_detect(contrast, "clr.Cardiomyocytes") ~ "Cardiomyocytes",
           str_detect(contrast, "clr.Fibroblast") ~ "Fibroblasts"
         )) |>
  arrange(desc(qscore.adj)) |> 
  slice_head(n = n.terms) |> 
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
    title = md("Cell-type adjusted GO terms associated with WT and cmAKO responses to MI"),
  ) |>  
  opt_row_striping() |> 
  cols_width(Description ~ 240) |> 
  tab_style(
    style = list(
      align = "center",
      cell_fill("grey"),
      cell_text(color = "black", weight = "bold")),
    locations = cells_row_groups())

file = "results/11_clusterProfiler/tables/supp_4a_top20_mi_interaction"
gtsave(tab, paste0(file, ".html"))
webshot::webshot(url = paste0(file, ".html"), file = "results/supp_figs/4d_top20_mi_interaction.png", 
                 vwidth = 700, vheight = n.terms*100, zoom = 3)



#### Supp 4e ####
n.terms <- 20
tab <- go |> 
  filter(qscore.adj >= 1.3 & contrast %in% c("clr.Cardiomyocytes", "clr.Fibroblast")) |> 
  group_by(contrast) |> 
  mutate(Description = stringr::str_to_title(Description),
         delta_q = qscore.adj- qscore.unadj,
         contrast = case_when(
           str_detect(contrast, "MI_vs_Sham") ~ "Myocardial Infarction",
           str_detect(contrast, "cmAKO_vs_WT") ~ "cmAKO",
           str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "cmAKO:MI",
           str_detect(contrast, "clr.Cardiomyocytes") ~ "Cardiomyocytes",
           str_detect(contrast, "clr.Fibroblast") ~ "Fibroblasts"
         )) |>
  arrange(desc(qscore.adj)) |> 
  slice_head(n = n.terms) |> 
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
    title = md("GO terms associated with cardiomyocytes and fibroblasts")
  ) |>  
  opt_row_striping() |> 
  cols_width(Description ~ 300) |> 
  tab_style(
    style = list(
      align = "center",
      cell_fill("grey"),
      cell_text(color = "black", weight = "bold")),
    locations = cells_row_groups())


file = "results/11_clusterProfiler/tables/supp_4b_top20_celltypes"
gtsave(tab, paste0(file, ".html"))
webshot::webshot(url = paste0(file, ".html"), file = "results/supp_figs/4e_top20_celltype.png", 
                 vwidth = 700, vheight = n.terms*100, zoom = 3)



#### Save a tabbed xlsx of all GO terms

go.adj <- readRDS("data/processed/pathway_genesets/goadjusted_005.RDS")
go.unadj <- readRDS("data/processed/pathway_genesets/gounadjusted_005.RDS")

names.adj <- names(go.adj) 
go.adj <- lapply(names(go.adj),
                 function(x){
                   data <- go.adj[[x]]@result
                   data <- data |> mutate(qscore = -log(p.adjust, base=10))
                   data$contrast <- x
                   return(data)})

names(go.adj) <- names.adj
names.unadj <- names(go.unadj)
go.unadj <- lapply(names(go.unadj),
                   function(x){
                     data <- go.unadj[[x]]@result
                     data <- data |> mutate(qscore = -log(p.adjust, base=10))
                     data$contrast <- x
                     return(data)}) 
names(go.unadj) <- names.unadj


# Generate an Excel file with each result as a separate tab

library(openxlsx)

# Create a new workbook
  wb <- createWorkbook()
  
  ## unadjusted results
  # Treatment MI vs Sham
  addWorksheet(wb, "MI_vs_Sham_unadjusted")
  writeData(wb, "MI_vs_Sham_unadjusted", go.unadj$treatment_MI_vs_Sham)
  
  # Gene x MI
  addWorksheet(wb, "GenotypeXTreatment_unadjusted")
  writeData(wb, "GenotypeXTreatment_unadjusted", go.unadj$treatmentMI.genotypecmAKO)
  
  # Genotype
  addWorksheet(wb, "Genotype_unadjusted")
  writeData(wb, "Genotype_unadjusted", go.unadj$genotype_cmAKO_vs_WT)
  
  ## Cell-type adjusted results
  # Treatment MI vs Sham
  addWorksheet(wb, "MI_vs_Sham_adjusted")
  writeData(wb, "MI_vs_Sham_adjusted", go.adj$treatment_MI_vs_Sham)
  
  # Gene x MI
  addWorksheet(wb, "GenotypeXTreatment_adjusted")
  writeData(wb, "GenotypeXTreatment_adjusted", go.adj$treatmentMI.genotypecmAKO)
  
  # Genotype
  addWorksheet(wb, "Genotype_adjusted")
  writeData(wb, "Genotype_adjusted", go.adj$genotype_cmAKO_vs_WT)
  
  # Cardiomyocytes
  addWorksheet(wb, "Cardiomyocytes_adjusted")
  writeData(wb, "Cardiomyocytes_adjusted", go.adj$clr.Cardiomyocytes)
  
  # Fibroblasts
  addWorksheet(wb, "Fibroblasts_adjusted")
  writeData(wb, "Fibroblasts_adjusted", go.adj$clr.Fibroblast)
  
  # Save the workbook
  saveWorkbook(wb, "results/supp_data/go_terms_all.xlsx", overwrite = TRUE)


