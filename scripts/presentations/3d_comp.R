library( rgl )
library(magick)

# Load libs
libs <- c("tidyverse", "RColorBrewer", "reshape2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#### Loading and formatting of data ####
# Load compositions
decon.whole <- read.csv("data/processed/compositions/whole_samples.csv")

#### Plot whole samples ####
#brewer.pal(n=8,"Paired")

my_palette <- c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00")

# Use the ordered new.id values to reorder the original data frame
decon.plotly <- decon.whole |> 
  mutate(Genotype_Treatment = paste(genotype, "-", treatment)) |> 
  mutate(Genotype_Treatment = factor(Genotype_Treatment, levels = c("WT - Sham", "WT - CAD", "cmAKO - Sham", "cmAKO - CAD"))) |> 
  dplyr::select(new.id, CellType, Prop, Genotype_Treatment) |> 
  pivot_wider(names_from = CellType, values_from = Prop)

fig <- plotly::plot_ly(decon.plotly, x = ~Cardiomyocytes, y= ~Fibroblast, 
                       z= ~Macrophage, color = ~Genotype_Treatment, colors = my_palette)

fig <- fig %>% plotly::add_markers()
fig <- fig %>% plotly::layout(
                scene = list(xaxis = list(title = 'Cardiomyocytes'),
                                           yaxis = list(title = 'Fibroblast'),
                                           zaxis = list(title = 'Macrophage'),
                                           camera = 
                                            list(eye = list(x = -1.25, y = 1.25, z = 1.25))))

fig


# save the widget
library(htmlwidgets)
saveWidget(fig, file="results/7_plot_comps/3d_comps.html")
