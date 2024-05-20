
# List libraries
libs <- c("Seurat", "SeuratDisk") # list libraries here
# Require all of them
lapply(libs, require, character.only = T)

rm(libs)

sn <- LoadH5Seurat("data/processed/single_cell/celltype_labeled.h5seurat")


# Function to make a data frame of cell type ratios
simulate_ratios <- function(major_cell, major_prop, 
                            cell_types, range, step_size, replicates, noise) {
  
  # get group tiers for major cell type
  # these are the percentages of the major cell type for which groups will be centered
  major_props_groups <- seq(major_prop - range, major_prop + range, step_size)
  
  # get values for major cell type
  # each value in the major_props_groups vector is used as a mean with n replicates pulled from a normal distribution with noise
  major_props <- lapply(major_props_groups, function(x){rnorm(replicates, x, noise)}) |>
    unlist()
  
  # Set major cell type and proportion
  minor <- cell_types[!(cell_types %in% c(major_cell))]
  
  
  # Define the range for proportions
  minor_props <- rep((1-major_props)/length(minor), length(minor)) |>
    matrix(ncol = length(minor)) |>
    as.data.frame()
  colnames(minor_props) <- minor
  
  # set up props dataframe
  combos <- c(major_props) |>
    matrix(ncol = 1) |>
    as.data.frame()
  colnames(combos) <- c(major_cell)
  
  # add minor cell types to props data frame
  props <- cbind(combos, minor_props)
  props$pct.change <- as.factor(rep(major_props_groups - major_prop, each = replicates))
  props$pct.change <- relevel(props$pct.change, ref = "0")
  rownames(props) <- paste0("mix_", 1:nrow(props))
  props[,2] <- (props[,1] + rowMeans(props[,3:length(cell_types)])) / 4
  props[,3] <- (props[,1] + rowMeans(props[,4:length(cell_types)])) / 5
  props[,2:length(cell_types)] <- sapply(2:length(cell_types), function(x){props[,x] * sample(rnorm(1000, 1, noise), length(major_props))})
  
  # Normalize back to proportional sum == 1
  props[,1:length(cell_types)]  <- props[,1:length(cell_types)]/rowSums(props[,1:length(cell_types)])
  
  return(props)
}

# Define arguments
major.cell <- "Cardiomyocytes"
major.prop <- 0.5
cell.types <- Idents(sn) |> unique() |> as.character()
range <- 0.3
step.size <- 0.1
replicates <- 5
noise <- 0.01

ratios <- simulate_ratios(major_cell = major.cell, 
                          major_prop = major.prop, 
                          cell_types = cell.types, 
                          range = range, 
                          step_size = step.size, 
                          replicates = replicates,
                          noise = noise)


# Function to aggregte the gene expression of cells from the same type into a single propability vector
getCellProfile <- function(sn, cell_type){
  # Get the cell.ids that are the cell type of interest
  cell.ids  <- Idents(sn)[which(Idents(sn) == cell_type)] |>
    names() 
  
  # sum expression of cells with those ids
  cell.type.profile <- rowSums(sn@assays$RNA@counts[,cell.ids]) |>
    as.numeric()  

  # Get the proportions of counts of each gene relative to the sum of all genes 
  cell.type.profile <- cell.type.profile/sum(cell.type.profile)
  assign(as.character(cell_type), cell.type.profile) # name the profile with the cell type
  return(get(as.character(cell_type)))
}

# make a df that has cell types as columns, genes as rows, and probablity of a single gene being expressed by a given cell type in each matrix cell
cell.profiles <- sapply(cell.types, function(x){getCellProfile(sn, x)}) |> as.data.frame()
row.names(cell.profiles) <- sn@assays$RNA@counts |> row.names() 

sampleExpression <- function(sample, umi_target){
  cell.umi <- ratios[sample, cell.types] * umi_target 
  # Multiply each column in cell.profiles by the corresponding value in cell.umi
  cell.prof <- cell.profiles * cell.umi[col(cell.profiles)]
  sample.exp <- rowSums(cell.prof) |> round()
  sample.exp
}

# Make a counts matrix
sim.counts <- sapply(row.names(ratios), function(x){sampleExpression(x, 25000000)})

# Save counts and ratios
write.csv(sim.counts, "data/processed/deseq_simulation/simulated_counts.csv")

write.csv(ratios, "data/processed/deseq_simulation/simulated_ratios.csv")