# Pathfinder
libs <- c("tidyverse", "pathfindR", "biomaRt", "DESeq2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Load the results and expression matrix
raw.res <- readRDS("data/processed/models/adjusted_de_interaction.RDS")
jensen.bulk <- read.csv("data/processed/bulk/jensen_bulk_clean.csv",  row.names = 1, check.names = F)

### Get mouse pathways ####
# Import them if they don't yet exist
if(!file.exists("data/processed/pathway_genesets/mmu_kegg_genes.RDS") |
   !file.exists("data/processed/pathway_genesets/mmu_kegg_descriptions.RDS")){
  gsets_list <- get_gene_sets_list(
    source = "KEGG",
    org_code = "mmu"
  )
  mmu_kegg_genes <- gsets_list$gene_sets
  mmu_kegg_descriptions <- gsets_list$descriptions

  saveRDS(mmu_kegg_genes, "data/processed/pathway_genesets/mmu_kegg_genes.RDS")
  saveRDS(mmu_kegg_descriptions, "data/processed/pathway_genesets/mmu_kegg_descriptions.RDS")
}else{
# Load them if they do exist
  mmu_kegg_genes <- readRDS("data/processed/pathway_genesets/mmu_kegg_genes.RDS")
  mmu_kegg_descriptions <- readRDS("data/processed/pathway_genesets/mmu_kegg_descriptions.RDS")
}

# Get STRING PIN
if(!file.exists("data/processed/pathway_genesets/mmusculusPIN.sif")){
## Downloading the STRING PIN file to tempdir
url <- "https://stringdb-downloads.org/download/protein.links.v12.0/10090.protein.links.v12.0.txt.gz"
path2file <- file.path(tempdir(check = TRUE), "STRING.txt.gz")
download.file(url, path2file)

## read STRING pin file
mmu_string_df <- read.table(path2file, header = TRUE)

## filter using combined_score cut-off value of 800
mmu_string_df <- mmu_string_df[mmu_string_df$combined_score >= 800, ]

## fix ids
mmu_string_pin <- data.frame(
  Interactor_A = sub("^10090\\.", "", mmu_string_df$protein1),
  Interactor_B = sub("^10090\\.", "", mmu_string_df$protein2)
)


# library(biomaRt)

mmu_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

converted <- getBM(
  attributes = c("ensembl_peptide_id", "mgi_symbol"),
  filters = "ensembl_peptide_id",
  values = unique(unlist(mmu_string_pin)),
  mart = mmu_ensembl
)
mmu_string_pin$Interactor_A <- converted$mgi_symbol[match(mmu_string_pin$Interactor_A, converted$ensembl_peptide_id)]
mmu_string_pin$Interactor_B <- converted$mgi_symbol[match(mmu_string_pin$Interactor_B, converted$ensembl_peptide_id)]
mmu_string_pin <- mmu_string_pin[!is.na(mmu_string_pin$Interactor_A) & !is.na(mmu_string_pin$Interactor_B), ]
mmu_string_pin <- mmu_string_pin[mmu_string_pin$Interactor_A != "" & mmu_string_pin$Interactor_B != "", ]

# remove self interactions
self_intr_cond <- mmu_string_pin$Interactor_A == mmu_string_pin$Interactor_B
mmu_string_pin <- mmu_string_pin[!self_intr_cond, ]

# remove duplicated inteactions (including symmetric ones)
mmu_string_pin <- unique(t(apply(mmu_string_pin, 1, sort))) # this will return a matrix object

mmu_string_pin <- data.frame(
  A = mmu_string_pin[, 1],
  pp = "pp",
  B = mmu_string_pin[, 2]
)

path2SIF <- file.path("data/processed/pathway_genesets", "mmusculusPIN.sif")
write.table(mmu_string_pin,
            file = path2SIF,
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE
)
}
path2SIF <- file.path("data/processed/pathway_genesets", "mmusculusPIN.sif")

# Run pathfindR ####
# Pull out interaction term results
res.raw <- results(raw.res, name="treatmentTAC.genotypeKO") |> as.data.frame()
res.raw <- data.frame(genes = row.names(res.raw),
                      log2FC = res.raw$log2FoldChange,
                      padj = res.raw$padj) |> 
  subset(!is.na(padj))
 
example_mmu_output <- run_pathfindR(
  input = res.raw,
  convert2alias = FALSE,
  gene_sets = "Custom",
  custom_genes = mmu_kegg_genes,
  custom_descriptions = mmu_kegg_descriptions,
  pin_name_path = path2SIF,
  p_val_threshold = 0.01
)



enrichment_chart(example_mmu_output,  top_terms = 50, plot_by_cluster = T)
clustered_fuzzy <- cluster_enriched_terms(example_mmu_output, method = "fuzzy")

enrichment_chart(clustered_fuzzy)

example_mmu_output
term_gene_heatmap(result_df = example_mmu_output, genes_df = res.raw)

UpSet_plot(result_df = example_mmu_output, genes_df = res.raw)
