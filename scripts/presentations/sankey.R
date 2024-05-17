# Sankey Plot

# Load libs
libs <- c("tidyverse", "clusterProfiler","patchwork","data.table") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

# Check if the simplified GO output exists, make it if not
if(!file.exists("data/processed/pathway_genesets/go_all_simp_005.RDS")){
# Load data
go.adj <- readRDS("data/processed/pathway_genesets/goadjusted_005.RDS")
go.unadj <- readRDS("data/processed/pathway_genesets/gounadjusted_005.RDS")

# Make names unique and join into one obj
names(go.adj) <- paste(names(go.adj), "adj", sep = "_")
names(go.unadj) <- paste(names(go.unadj), "unadj", sep = "_")
go <- c(go.adj, go.unadj)
rm(go.adj, go.unadj)

# Simplify terms by semantic similarity
go.simp <- lapply(go, clusterProfiler::simplify, cutoff = 0.6)
saveRDS(go.simp, "data/processed/pathway_genesets/go_all_simp_005.RDS")
}else{
go.simp <- readRDS("data/processed/pathway_genesets/go_all_simp_005.RDS")
}

#### Sankey Diagram ####

# Library
library(networkD3)
library(dplyr)

# Merge go datasets
# Make into dataframes of result with column for contrast and model version
go.merged <- lapply(names(go.simp),
                       function(x){
                         data <- go.simp[[x]]@result
                         data$contrast <- x
                         return(data)}) |> 
                    rbindlist() |> 
                    as.data.frame() |> 
                    select(ID, Description, p.adjust, contrast)

# Make a before and after column from adj/unadj
go.merged <- go.merged |> 
  mutate(
    model = case_when(
      str_detect(contrast, "_unadj") ~ "unadjust",
      .default = "adjusted"
    ),
    contrast = case_when(
      str_detect(contrast, "treatment_MI_vs_Sham") ~ "MI",
      str_detect(contrast, "genotype_cmAKO_vs_WT") ~ "cmAKO",
      str_detect(contrast, "treatmentMI.genotypecmAKO") ~ "MI:cmAKO",
      str_detect(contrast, "Fibroblast") ~ "Fibroblast",
      str_detect(contrast, "Cardio") ~ "Cardiomyocytes"
    )
  )


combinations <- expand.grid(unique(go.merged$model), unique(go.merged$contrast))
combinations <- combinations[c(1:5,7,9,10),] # hardcoded to leave out variables not used before adjustment in the preadjustmnet combinations
combination.list <-c()
# Make a loop to pull the IDs from each combination
for(i in 1:nrow(combinations)){
  combination.list[[i]] <- go.merged |> filter(model ==  combinations[i,1] & contrast == combinations[i,2]) |> pull(ID) 
}

names(combination.list) <- paste(combinations$Var1, combinations$Var2, sep = "_")

#combination.list[["bait"]] <- unlist(combination.list)
lengths <- lapply(names(combination.list), function(x){length(combination.list[[x]])})
lengths <- data.frame(source = names(combination.list),target = names(combination.list), length = unlist(lengths))

# Make the combination table
comb.table <- crossprod(table(stack(combination.list))) |> as.data.frame()
comb.table$source <- row.names(comb.table)

comb.table <-pivot_longer(comb.table, cols = names(combination.list), values_to = "value", names_to = "target")

# add the variable totals, make lost variable in target and gained in source
#comb.table <- left_join(comb.table, lengths)
comb.table |> 
  group_by(source) |> 
  mutate(lost = max(value) - (sum(value) - max(value))) |> 
  subset(!is.na(length))


comb.table <- comb.table |> mutate(
  filt = case_when(
    str_detect(source, "unadjust") & str_detect(target, "adjusted") ~ TRUE,
    .default = FALSE),
  filt = case_when(
    target == source ~ TRUE,
    .default = filt)
    ) 

comb.table <- comb.table |> filter(filt == TRUE) 

comb.table <- comb.table |> mutate(
  source = case_when(
    target == source & str_detect(source, "adjusted") ~ "Gained after adjusting",
    .default = source),
  target = case_when(
    target == source & str_detect(target, "unadjust") ~ "Lost once Adjusted",
    .default = target))

# Want to add a middle row of nodes so that Gained and Lost nodes start/end there
# subset to the non gained/lost rows
comb.target <- comb.table |> 
  mutate(target = str_replace(target, ".*_", "")) |> 
  filter(source != "Gained after adjusting" & target != "Lost once Adjusted") |> 
  mutate(target = paste(target, "middle", sep = "_"))

comb.source <- comb.table |> 
  mutate(source = str_replace(source, ".*_", "")) |> 
  filter(source != "Gained after adjusting" & target != "Lost once Adjusted") |> 
  mutate(source = paste(source, "middle", sep = "_"))

comb.change <- comb.table |> 
  filter(source == "Gained after adjusting" | target == "Lost once Adjusted")

comb.all <- rbind(comb.change, comb.target) |> rbind(comb.source)
# duplicate them - for one, change target to be target_middle, for the other, make source into source_middle
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(comb.all$source), 
         as.character(comb.all$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
comb.all$IDsource <- match(comb.all$source, nodes$name)-1 
comb.all$IDtarget <- match(comb.all$target, nodes$name)-1

# Make the Network
p <- sankeyNetwork(Links = comb.all, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)
p


####################
#### Rachel's sankey #####
########################
library(networkD3)

# List the filter steps
nodes <- c("WGS Demo ","EHR",
           "EU Ancestry","PSS Survey",
           "Controls","Cases")
nodes <- c("All samples", nodes) # Add a node to indicate that the first number is all of the samples

# group sizes
values <- c(245388,184536,100504,62070,27304,24144,3160)

# Don't include the "All Samples" # in the links 
links <- data.frame("value" = values[-1])

# This was manual. For each node, where does it start and where does that value go?
# Excluding "All Samples", start from 0
links$source <- c(0, 1, 2, 3, 4, 4)
links$target <- c(1,2,3,4,5,6)

# This was my solution to adding lost samples
# its the # of steps with losses
n.delta <- 4
# Get lost samples values
deltas <- data.frame(value = rep(NA, n.delta),
                     source = rep(NA, n.delta),
                     target = rep(NA, n.delta))

# Make a new df to make the lost values in
j <- length(nodes)
for(i in 1:n.delta){
  deltas[i, "value"] <- values[i] - values[i+1]
  deltas[i, "source"] <- i-1
  deltas[i, "target"] <- j
  j <- j + 1
}  

# Add new node names for the lost values
nodes <- c(nodes, paste("Lost in", nodes[2:(n.delta+1)], "Step"))

# Join the original and new dfs
links <- rbind(links, deltas)

# Add numbers to the labels
# the "\n" new line delimiter didn't work :\
nodes <- paste0(nodes, "\n(n = ", links$value, ")")

nodes <- data.frame(name = nodes)

sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize= 20, nodeWidth = 40, sinksRight = F,
              colourScale =JS("d3.scaleOrdinal(d3.schemeCategory20);"))

