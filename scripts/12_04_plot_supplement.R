# This script is meant to generate a figure which shows the performance of 
# different representations of cell type composition in simulated differential
# gene expression analysis

## Retrieve data ####

# List libraries
libs <- c("readr", "purrr", "dplyr", "stringr", "ggplot2", "grid") # list libraries here
lapply(libs, require, character.only = T) # Require all of them
rm(libs)

# Load data
file.path <- "data/processed/deseq_simulation/batched_output/"
files <- paste0(file.path, list.files(file.path))
all_data <- files %>%
  map_df(~read_csv(., show_col_types = FALSE)) # read and join csv files

# Format major cell type percents as numbers
all_data$pct.change <- all_data$pct.change |> str_replace( "_", "-") |> as.numeric()

# Make a palette for the figure
pal <- c("#969696",
         "#c51b8a",
         "#006837",
         "#fd8d3c"
)

# Modify legend text 

models.legend <- list(
  "1"    = ~ "No cell types",
  "2"    = ~ "Untransformed CM proportion",
  "3"    = ~ "CLR-transformed CM proportion",
  "4"    = ~ "PC1 of cell type proportions"
)


## Plot the figure ####

p.sim <- all_data |> 
  mutate(pct.sig = pct.sig*100,
         pct.change = 50 + 100*pct.change,
         model = case_when(
           str_detect(model, "1") ~  "Cell-type blind",
           str_detect(model, "2") ~  "Untransformed CM proportion",
           str_detect(model, "3") ~  "CLR-transformed CM proportion",
           str_detect(model, "4") ~  "PC1 of cell type proportions"
         )) |>
  #filter(model %in% c(2:4)) |> 
  ggplot(aes(x = pct.change, y = pct.sig, color = as.factor(model))) +
  geom_point(alpha = 0.2, size = 0.6) +
  geom_smooth(se = F, method = "loess", size = 1, alpha = 0.2, span = 0.3) +
  facet_wrap(~model, nrow = 2) +
  scale_color_manual(values = pal,
                     name = "DESeq2 model variables",
                     #labels = models.legend,
                     #breaks = names(models.legend)
                     ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 8),
    legend.text = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
    legend.key.size = unit(0.5, "lines"),  # Adjust size of the keys
    legend.spacing.x = unit(0.5, "mm"),    # Reduce spacing between columns of the legend
    legend.spacing.y = unit(0.5, "mm"),    # Reduce spacing between rows of the legend
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_line(color = "grey"),
    plot.margin = unit(c(0,0,0,0), units = "cm"),# remove minor gridlines
  ) +
  coord_cartesian(clip = 'off') +
  labs(y= "% DEGs vs. 50% CMs", x = "Cardiomyocyte Percentage of Samples") +
  guides(color = guide_legend(nrow = 4))
#Maybe add cm% bar under plot
# make legend text more direct/clear

#Save outputs
if(!dir.exists("results/supp_figs")){
  dir.create("results/supp_figs")
}

png(file = "results/supp_figs/4a_sim_de.png",
    width = 3.0603*1.5, 
    height = 1.9348* 1.5,
    units = "in",
    res = 600)

p.sim

dev.off()

