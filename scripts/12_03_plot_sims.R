# List libraries
libs <- c("readr", "purrr", "dplyr", "stringr", "ggplot2", "grid", "tidyr", "tidytext") # list libraries here
# Require all of them
lapply(libs, require, character.only = T)

rm(libs)

files <- list.files("data/processed/deseq_simulation/batched_output/")

# temp subset to 1:2_1:2
#temp.df <- expand_grid(as.character(3), paste0("_", 1:1))
#temp.files <- paste0(temp.df[[1]], temp.df[[2]], ".csv")
files <- paste0("data/processed/deseq_simulation/batched_output/", files)


all_data <- files %>%
  map_df(~read_csv(., show_col_types = FALSE))

pct_models <- all_data |> 
  filter(str_detect(pct.change, "0.")) |> 
  select(-pct.sig)

pct_models$pct.change <- pct_models$pct.change |> str_replace( "_", "-") |> as.numeric()


pal <- c("#969696",
         "#a6cee3",
         "#1f78b4",
         "#fdbf6f",
         "#ff7f00",
         "#33a02c"
)


# Make list of models
models <- list(
  unadjusted    = ~ 0 + pct.change,
  raw_cm        = ~ 0 + pct.change + Cardiomyocytes,
  raw_cm_fb     = ~ 0 + pct.change + Cardiomyocytes + Fibroblast,
  clr_cm        = ~ 0 + pct.change + clr_Cardiomyocytes,
  clr_cm_fb     = ~ 0 + pct.change +  clr_Cardiomyocytes + clr_Fibroblast,
  pc1           = ~ 0 + pct.change + PC1
)
# Modify legend text 
models.legend <- list(
  "1"    = ~ "No cell types",
  "2"    = ~ "Raw prop 1x",
  "3"    = ~ "Raw prop 2x",
  "4"    = ~ "CLR 1x",
  "5"    = ~ "CLR 2x",
  "6"    = ~ "PC1"
)


# Calculate metrics and add them to the tibble.
metrics_df <- pct_models %>%
  mutate(
    precision = `True positive` / (`True positive` + `False positive`),
    recall = `True positive` / (`True positive` + `False negative`),
    f1.score = 2 * (precision * recall) / (precision + recall),
    accuracy = (`True positive` + `True negative`) /
      (`True positive` + `False positive` + `False negative` + `True negative`),
    specificity = `True negative` / (`True negative` + `False positive`)
  )


pct_long <- pct_models |> pivot_longer(cols = c("True negative",  "False positive", "False negative","True positive"))

p.sim <- metrics_df |> 
  mutate(pct.change = 50 + 100*pct.change) |> 
  ggplot(aes(x = pct.change, y = f1_score, color = as.factor(model))) +
  geom_point(alpha = 0.2, size = 0.6) +
  geom_smooth(se = F, method = "loess", linewidth = 1, alpha = 0.2, span = 0.3) +
  scale_color_manual(values = pal,
                     name = "",
                     labels = models.legend,
                     breaks = names(models.legend)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    text = element_text(size = 8),
    legend.text = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
    legend.key.size = unit(0.5, "lines"),  # Adjust size of the keys
    legend.spacing.x = unit(0.5, "mm"),    # Reduce spacing between columns of the legend
    legend.spacing.y = unit(0.5, "mm"),    # Reduce spacing between rows of the legend
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_line(color = "grey"),
    plot.margin = unit(c(0.2,0,0,0), units = "cm"),
    rect = element_rect(fill = "transparent")# remove minor gridlines
  ) +
  coord_cartesian(clip = 'off') +
  labs(y= "F1 Score (higher is more accurate)", x = "Simulated Percent Cardiomyocytes in Samples",
       subtitle = str_wrap("F1 scores balance recall (% DEGs not due to composition) and precision (% true DEGs the model detected)",
                          width = 80)) +
  guides(color = guide_legend(nrow = 1))
#Maybe add cm% bar under plot
# make legend text more direct/clear

#Save outputs
if(!dir.exists("results/12_sim_de")){
  dir.create("results/12_sim_de")
}

png(file = "results/12_sim_de/plot_sims.png",
    width = 5, 
    height = 3 ,
    units = "in",
    res = 600)

p.sim

dev.off()



# Calculate metrics and add them to the tibble.
metrics_df <- pct_models %>%
  mutate(
    precision = `True positive` / (`True positive` + `False positive`),
    recall = `True positive` / (`True positive` + `False negative`),
    f1.score = 2 * (precision * recall) / (precision + recall),
    accuracy = (`True positive` + `True negative`) /
      (`True positive` + `False positive` + `False negative` + `True negative`),
    specificity = `True negative` / (`True negative` + `False positive`)
  )


metrics_trimmed <- metrics_df |> 
  filter(model %in% c("1", "2", "4", "6"))

models.legend <- list(
  "1"    = "No cell types",
  "2"    = "Raw prop",
  "3"    = "Raw prop 2x",
  "4"    = "CLR",
  "5"    = "CLR 2x",
  "6"    = "PC1"
)


for(i in names(models.legend)){
metrics_df$model[which(as.character(metrics_df$model) == i)] <- models.legend[[i]]
}

# check how many iterations there are for each group/bin
temp <- metrics_df |> 
  mutate(bin = case_when(
    abs(pct.change) >= 0.1 ~ "Major change",
    abs(pct.change) <  0.1 ~ "Minor change"))

table(temp$model, temp$bin)

# Make plotable sumarry table and plot
summaryTable <- metrics_df |> 
  mutate(bin = case_when(
    abs(pct.change) >= 0.1 ~ "Major change (>10% shift)",
    abs(pct.change) <  0.1 ~ "Minor change (<10% shift)")) |> 
  group_by(bin, model) %>%
  summarise(
    across(precision:specificity,
           list(mean = ~mean(., na.rm = TRUE),
                stdev = ~sd(., na.rm = TRUE)),
           .names = "{.col}_{.fn}")
  ) %>%
  pivot_longer(
    cols = -c(bin, model),
    names_to = c("metric", ".value"),
    names_sep = "_"
  ) |> 
  filter(metric %in% "f1.score")

f1Plot <- ggplot(summaryTable, aes(x = reorder_within(model, mean, bin), y = mean, color = model)) +
  geom_errorbar(aes(ymin = mean - stdev, ymax = mean + stdev),
                width = 0.2, size = 1) +
  geom_point(size = 3) +
  facet_grid(~bin, scales = "free_x") +
  scale_x_reordered() +
  labs(x = "Cell type representation in model",
       y = "F1 Score",
       color = "Model") +
  scale_color_manual(values = c("#66c2a5", "grey", "#fc8d62", "#e78ac3")) +
  theme_minimal() + 
  theme(plot.title = element_blank(),
        legend.position = "none")


png(file = "results/12_sim_de/f1_plot.png",
    width = 5, 
    height = 3 ,
    units = "in",
    res = 600)

f1Plot

dev.off()


# plot results associated with the cell type variable

# rename the columns to reflect that false postives are actually true postives since we're looking at cell-type assocaited changes


cell_models <- all_data |> 
  filter(!str_detect(pct.change, "0."))

colnames(cell_models) <- colnames(cell_models)[c(3, 4, 1, 2, 5:9)]


cell_metrics <- cell_models %>%
  mutate(
    precision = `True positive` / (`True positive` + `False positive`),
    recall = `True positive` / (`True positive` + `False negative`),
    f1_score = 2 * (precision * recall) / (precision + recall),
    accuracy = (`True positive` + `True negative`) /
      (`True positive` + `False positive` + `False negative` + `True negative`),
    specificity = `True negative` / (`True negative` + `False positive`)
  )


data_long <- cell_models |> pivot_longer(cols = c("True negative",  "False positive", "False negative","True positive"))

p.sim <- cell_metrics |>  
  ggplot(aes(x = model, y = f1_score, color = as.factor(model))) +
  geom_point(alpha = 0.2, size = 0.6) +
  scale_color_manual(values = pal,
                     name = "",
                     labels = models.legend,
                     breaks = names(models.legend)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    text = element_text(size = 8),
    legend.text = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
    legend.key.size = unit(0.5, "lines"),  # Adjust size of the keys
    legend.spacing.x = unit(0.5, "mm"),    # Reduce spacing between columns of the legend
    legend.spacing.y = unit(0.5, "mm"),    # Reduce spacing between rows of the legend
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_line(color = "grey"),
    plot.margin = unit(c(0,0,0,0), units = "cm"),
    rect = element_rect(fill = "transparent")# remove minor gridlines
  ) +
  coord_cartesian(clip = 'off') +
  labs(y= "% DEGs vs. 50% CMs", x = "Simulated Percent Cardiomyocytes in Samples") +
  guides(color = guide_legend(nrow = 1))
#Maybe add cm% bar under plot
# make legend text more direct/clear

#Save outputs
if(!dir.exists("results/12_sim_de")){
  dir.create("results/12_sim_de")
}

png(file = "results/12_sim_de/plot_sims_cellCovariates.png",
    width = 5, 
    height = 3 ,
    units = "in",
    res = 600)

p.sim

dev.off()
