# This script is meant to save the full DESeq2 results as a supplemental file
# There are two RDS files (one for each DESeq2 model: cell type aware and blind)


# Visualize DE 
libs <- c("tidyverse", "DESeq2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load the results and expression matrices
res.adj <- readRDS("data/processed/models/adjusted_de_interaction.RDS")
names <- resultsNames(res.adj)[-1]

# Store each results contrast as df in a list
res.adj <- lapply(names, function(x){
  results(res.adj, name=x) |> 
    as.data.frame() |> 
    rownames_to_column(var = "gene") 
})
names(res.adj) <- paste0(names, ".adj")

# Repeat for unadjusted results, then combine them
res.unadj <- readRDS("data/processed/models/unadjusted_de_interaction.RDS")
names <- resultsNames(res.unadj)[-1]

# Store each results contrast as df in a list
res.unadj <- lapply(names, function(x){
  results(res.unadj, name=x) |> 
    as.data.frame() |> 
    rownames_to_column(var = "gene")
})

names(res.unadj) <- paste0(names, ".unadj")

# Combine
res <- c(res.adj, res.unadj)

rm(res.unadj)
rm(res.adj)


## Generate CSVs for lost and gained sig genes ####
# > 10x change in padj
# Just the gene x treat interaction and MI contrasts

# Helper function to determine significant gain/reduction
find_significant_changes <- function(df_adj, df_unadj) {
  df <- inner_join(df_adj, df_unadj, by = "gene", suffix = c(".adj", ".unadj"))
  
  df <- df %>%
    mutate(padj_ratio = padj.unadj / padj.adj,
           significance_change = case_when(
             padj_ratio > 10 ~ "gain",      # unadjusted padj is 10x greater than adjusted
             padj_ratio < 0.1 ~ "reduce",   # adjusted padj is 10x greater than unadjusted
             TRUE ~ "no_change"),
           neg_log_10_padj.adj = -log10(padj.adj),
           neg_log_10_padj.unadj = -log10(padj.unadj))
  
  # Filter for gain and reduce significance changes
  gain <- df %>% filter(significance_change == "gain")
  reduce <- df %>% filter(significance_change == "reduce")
  
  list(gain = gain, reduce = reduce)
}

# Make data directory to store these CSVs
data.path <- "results/supp_data/DEGs"

if(!dir.exists(data.path)){
  dir.create(data.path)
}

# Function to generate CSV files
generate_csv_files <- function(res) {
  # Treatment MI vs Sham
  result_treatment <- find_significant_changes(res$treatment_MI_vs_Sham.adj, res$treatment_MI_vs_Sham.unadj)
  write.csv(result_treatment$gain, paste(data.path,"DEGs_treat_gain.csv", sep = "/"), row.names = FALSE)
  write.csv(result_treatment$reduce, paste(data.path,"DEGs_treat_reduce.csv", sep = "/"), row.names = FALSE)
  
  # Treatment MI genotype cmAKO
  result_genotype <- find_significant_changes(res$treatmentMI.genotypecmAKO.adj, res$treatmentMI.genotypecmAKO.unadj)
  write.csv(result_genotype$gain, paste(data.path, "DEGs_gene_treat_gain.csv", sep = "/"), row.names = FALSE)
  write.csv(result_genotype$reduce, paste(data.path, "DEGs_gene_treat_reduce.csv", sep = "/"), row.names = FALSE)
}

# Generate the CSV files
generate_csv_files(res)

### Make a tabbed excel sheet

library(openxlsx)

# Function to generate an Excel file with each result as a separate tab
generate_excel_file <- function(res, file_path) {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Treatment MI vs Sham
  result_treatment <- find_significant_changes(res$treatment_MI_vs_Sham.adj, res$treatment_MI_vs_Sham.unadj)
  addWorksheet(wb, "MI_vs_Sham_Gain")
  writeData(wb, "MI_vs_Sham_Gain", result_treatment$gain)
  
  addWorksheet(wb, "MI_vs_Sham_Reduce")
  writeData(wb, "MI_vs_Sham_Reduce", result_treatment$reduce)
  
  # Treatment MI genotype cmAKO
  result_genotype <- find_significant_changes(res$treatmentMI.genotypecmAKO.adj, res$treatmentMI.genotypecmAKO.unadj)
  addWorksheet(wb, "GenotypeXTreatment_Gain")
  writeData(wb, "GenotypeXTreatment_Gain", result_genotype$gain)
  
  addWorksheet(wb, "GenotypeXTreatment_Reduce")
  writeData(wb, "GenotypeXTreatment_Reduce", result_genotype$reduce)
  
  # Save the workbook
  saveWorkbook(wb, file_path, overwrite = TRUE)
}

# Set the path for the Excel file
excel_file_path <-  "results/supp_data/DEGs_Results_most_changed.xlsx"

# Execute the function to generate the Excel file
generate_excel_file(res, excel_file_path)



#### Make tabbed xlsx for ALL DEGs ####
# Each tab is a contrast from each model iteration (cell type adjusted or not)

# Function to generate an Excel file with each result as a separate tab
generate_excel_file <- function(res, file_path) {
  # Create a new workbook
  wb <- createWorkbook()
  
  ## unadjusted results
  # Treatment MI vs Sham
  addWorksheet(wb, "MI_vs_Sham_unadjusted")
  writeData(wb, "MI_vs_Sham_unadjusted", res$treatment_MI_vs_Sham.unadj)
  
  # Gene x MI
  addWorksheet(wb, "GenotypeXTreatment_unadjusted")
  writeData(wb, "GenotypeXTreatment_unadjusted", res$treatmentMI.genotypecmAKO.unadj)
  
  # Genotype
  addWorksheet(wb, "Genotype_unadjusted")
  writeData(wb, "Genotype_unadjusted", res$genotype_cmAKO_vs_WT.unadj)
  
  ## Cell-type adjusted results
  # Treatment MI vs Sham
  addWorksheet(wb, "MI_vs_Sham_adjusted")
  writeData(wb, "MI_vs_Sham_adjusted", res$treatment_MI_vs_Sham.adj)
  
  # Gene x MI
  addWorksheet(wb, "GenotypeXTreatment_adjusted")
  writeData(wb, "GenotypeXTreatment_adjusted", res$treatmentMI.genotypecmAKO.adj)
  
  # Genotype
  addWorksheet(wb, "Genotype_adjusted")
  writeData(wb, "Genotype_adjusted", res$genotype_cmAKO_vs_WT.adj)
  
  # Cardiomyocytes
  addWorksheet(wb, "Cardiomyocytes_adjusted")
  writeData(wb, "Cardiomyocytes_adjusted", res$clr.Cardiomyocytes.adj)
  
  # Fibroblasts
  addWorksheet(wb, "Fibroblasts_adjusted")
  writeData(wb, "Fibroblasts_adjusted", res$clr.Fibroblast.adj)
  
  # Save the workbook
  saveWorkbook(wb, file_path, overwrite = TRUE)
}

# Set the path for the Excel file
excel_file_path <-  "results/supp_data/DEGs_Results_All.xlsx"

# Execute the function to generate the Excel file
generate_excel_file(res, excel_file_path)
