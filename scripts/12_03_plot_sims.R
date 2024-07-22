# List libraries
libs <- c("readr", "purrr", "dplyr", "stringr", "ggplot2", "grid") # list libraries here
# Require all of them
lapply(libs, require, character.only = T)

rm(libs)

files <- list.files("data/processed/deseq_simulation/batched_output/")

files <- paste0("data/processed/deseq_simulation/batched_output/", files)

all_data <- files %>%
  map_df(~read_csv(., show_col_types = FALSE))

all_data$pct.change <- all_data$pct.change |> str_replace( "_", "-") |> as.numeric()


pal <- c("#969696",
         "#c51b8a",
         "#006837",
         "#fd8d3c"
)

# Modify legend text 
models.legend <- list(
  "1"    = ~ "No cell types",
  "2"    = ~ "CM proportion",
  "3"    = ~ "Cell types in model",
  "4"    = ~ "PC1"
)


p.sim <- all_data |> 
  mutate(pct.sig = pct.sig*100,
         pct.change = 50 + 100*pct.change) |>
  filter(model %in% c(1,3)) |> 
  ggplot(aes(x = pct.change, y = pct.sig, color = as.factor(model))) +
  geom_point(alpha = 0.2, size = 0.6) +
  geom_smooth(se = F, method = "loess", size = 1, alpha = 0.2, span = 0.3) +
  scale_color_manual(values = pal,
                     name = "DESeq2 model variables",
                     labels = models.legend,
                     breaks = names(models.legend)) +
  theme_minimal() +
  theme(
    legend.position = c(0.5, 0.75),
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
  labs(y= "% DEGs vs. to 50% CMs", x = "Cardiomyocyte Percentage of Samples") +
  guides(color = guide_legend(nrow = 4))
#Maybe add cm% bar under plot
# make legend text more direct/clear

#Save outputs
if(!dir.exists("results/12_sim_de")){
  dir.create("results/12_sim_de")
}

png(file = "results/12_sim_de/plot_sims.png",
    width = 3.0603, 
    height = 1.9348 ,
    units = "in",
    res = 600)

p.sim

dev.off()

