# Make a phenotype df for bulk RNAseq cell type fractions

# list files
files <- list.files("data/raw/fastq")[str_detect(list.files("data/raw/fastq"), "B6_")]

# Split name details
files.split <- str_split(files, "_")

# Compile in df
files.df <- data.frame(colname = files)
files.df$sex <- lapply(files.split, "[[", 2) |> unlist()
files.df$cell.type <- lapply(files.split, "[[", 3) |> unlist()
files.df$replicate <- lapply(files.split, "[[", 5) |> unlist()

# Save
write.csv(files.df, "data/raw/rau_fractions/celltype_pheno.csv", row.names = F)
