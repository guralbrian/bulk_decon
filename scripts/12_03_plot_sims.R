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
         "#fa9fb5", "#c51b8a",
         "#78c679", "#006837",
         "#fd8d3c","#e6550d"
)

# Modify legend text 
models.legend <- list(
  "1"    = ~ "No cell types in model",
  "2"    = ~ "CM proportion included in model",
  "3"    = ~ "Group + CM prop + 1x prop of minor cells",
  "4"    = ~ "CLR transofr CM proportion in cell mixture",
  "5"    = ~ "Group + CLR of CMs + 1x CLR of minor cells",
  "6"    = ~ "Group + PC1",
  "7"    = ~ "Group + PC1 + PC2"
)


p.sim <- all_data |> 
  mutate(pct.sig = pct.sig*100,
         pct.change = 50 + 100*pct.change) |>
  group_by(pct.change, model ) |> 
  mutate(upper = max(pct.sig),
         lower = min(pct.sig)) |> 
  filter(model %in% c(1,2)) |> 
  ggplot(aes(x = pct.change, y = pct.sig, color = as.factor(model))) +
  geom_point(alpha = 0.2, size = 5) +
  geom_smooth(se = F, method = "loess", size = 1.5, alpha = 0.2, span = 0.3) +
  #geom_ribbon(aes(ymin = lower, ymax = upper, color = as.factor(model)), alpha = 0.5) +
  scale_color_manual(values = pal,
                     name = "   DESeq2\nmodel design",
                     labels = models.legend,
                     breaks = names(models.legend)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    text = element_text(size = 30),
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_line(color = "grey")  # remove minor gridlines
  ) +
  #annotate("segment", x =   40 , xend = 32, y = 10, yend = 10,
  #         colour = "black", arrow = arrow()) + annotate("text", x = 36, y = 13, label = "Less Cardiomyocytes") +
  #annotate("segment", x =   60 , xend = 68, y = 10, yend = 10,
  #         colour = "black", arrow = arrow()) + annotate("text", x = 64, y = 13, label = "More Cardiomyocytes") +
  coord_cartesian(clip = 'off') +
  labs(y= "% DEGs compared to 50% CMs group", x = "Cardiomyocyte Percentage of Samples")

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
