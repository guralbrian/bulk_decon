# The purpose of this script is to merge the indiviual GO results into a tabbed 
# xlsx sheet

## Load packages
library(openxlsx)
library(tidyverse)

## Load data

# Get file names, remove DEG part
files <- list.files("results/supp_data/go_terms/")

files <- files[str_detect(files, "DEGs")] # Take only the DEG files

res <- lapply(files, function(x){read.csv(paste0("results/supp_data/go_terms/",x))[,c(2:11)]}) # Load them in to a list

tab_names <- lapply(str_split(basename(files), "DEGs_"), "[[", 2) # Get file name

tab_names <- lapply(str_split(tab_names, ".csv"), "[[", 1) |> unlist() # Get rid of csv extension

names(res) <- tab_names


# Function to generate an Excel file with each result as a separate tab
generate_excel_file <- function(res, file_path) {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Treatment MI vs Sham
  addWorksheet(wb, "treat_gain")
  writeData(wb, "treat_gain", res$treat_gain)
  
  addWorksheet(wb, "treat_reduce")
  writeData(wb, "treat_reduce", res$treat_reduce)
  
  # Treatment MI x genotype cmAKO
  addWorksheet(wb, "gene_treat_gain")
  writeData(wb, "gene_treat_gain", res$gene_treat_gain)
  
  addWorksheet(wb, "gene_treat_reduce")
  writeData(wb, "gene_treat_reduce", res$gene_treat_reduce)
  
  # Save the workbook
  saveWorkbook(wb, file_path, overwrite = TRUE)
}

# Set the path for the Excel file
excel_file_path <-  "results/supp_data/go_terms_most_changed_DEGs.xlsx"

# Execute the function to generate the Excel file
generate_excel_file(res, excel_file_path)
