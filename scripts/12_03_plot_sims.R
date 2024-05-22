# List libraries
libs <- c("readr", "purrr", "dplyr", "stringr", "ggplot2") # list libraries here
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
  "1"    = ~ "Group",
  "2"        = ~ "Group + CM prop",
  "3"     = ~ "Group + CM prop + 1x prop of minor cells",
  "4"        = ~ "Group + CLR of CMs",
  "5"     = ~ "Group + CLR of CMs + 1x CLR of minor cells",
  "6"           = ~ "Group + PC1",
  "7"          = ~ "Group + PC1 + PC2"
)


all_data |> 
  mutate(pct.sig = pct.sig*100,
         pct.change = 50 + 100*pct.change) |> 
  ggplot(aes(x = pct.change, y = pct.sig, color = as.factor(model))) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = F, method = "loess", size = 1.5, alpha = 0.2, span = 0.5) +
  scale_color_manual(values = pal,
                     name = "Model design",
                     labels = models.legend,
                     breaks = names(models.legend)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 12),
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank()  # remove minor gridlines
  ) +
  labs(y= "DE genes with p < 0.05 (percent)", x = "Simulated Cardiomyocyte Percent")

