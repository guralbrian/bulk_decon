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
  "1"    = ~ "No cell types in model",
  "2"    = ~ "Unchanged CM proportion",
  "4"    = ~ "CLR transform of CM proportion",
  "6"    = ~ "PC1 of cell proportions"
)


p.sim <- all_data |> 
  mutate(pct.sig = pct.sig*100,
         pct.change = 50 + 100*pct.change) |>
  filter(model %in% c(1,2,4,6)) |> 
  ggplot(aes(x = pct.change, y = pct.sig, color = as.factor(model))) +
  geom_point(alpha = 0.2, size = 5) +
  geom_smooth(se = F, method = "loess", size = 1.5, alpha = 0.2, span = 0.3) +
  scale_color_manual(values = pal,
                     name = "      DESeq2\nmodel variables",
                     labels = models.legend,
                     breaks = names(models.legend)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    text = element_text(size = 30),
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_line(color = "grey")  # remove minor gridlines
  ) +
  coord_cartesian(clip = 'off') +
  labs(y= "% DEGs compared to 50% CMs group", x = "Cardiomyocyte Percentage of Samples") +
  guides(color = guide_legend(nrow = 2))
#Maybe add cm% bar under plot
# make legend text more direct/clear

#Save outputs
if(!dir.exists("results/12_sim_de")){
  dir.create("results/12_sim_de")
}

png(file = "results/12_sim_de/plot_sims.png",
    width = 16, 
    height = 9,
    units = "in",
    res = 600)

p.sim

dev.off()
