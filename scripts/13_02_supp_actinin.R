# This script is meant to generate a figure which visualizes the area occupied 
# by actinin in IHC/RNAscope, which is supplement figure 5a

# Load data
data <- read.csv("data/raw/rnascope/alpha1_ar_image_quant.csv")

# Load libraries
library(tidyverse)

# Pull gene, sample type, and number from sample name
# Reusing script used to make main figure so some of this is unnecessary
data <- data |> 
  mutate(
    channel = factor(channel),
    gene.group = case_when(
      str_detect(sample_name, "Actinin") ~ "Actinin",  
      str_detect(sample_name, "CD45") ~ "CD45",  
      str_detect(sample_name, "DDR2") ~ "DDR2"),
    region = case_when(
      str_detect(region, " A") | str_detect(region, "remote") ~ "remote",
      str_detect(region, " B") | str_detect(region, "infarct") ~ "infarct",
      str_detect(region, " C") | str_detect(region, "border") ~ "border",
      .default = region)) |> 
  mutate(channel = case_when(area.of.spots == 396 ~ "3",
                             .default = channel)) |>  # Fix typo from data input
  mutate(gene = case_when(
    gene.group == "Actinin" & channel == 2 ~ "Actinin",
    gene.group == "Actinin" & channel == 3 ~ "Pik3r1",
    gene.group == "Actinin" & channel == 4 ~ "Zbtb16"),
    number = number.of.spots.detected,
    area = area.of.spots,
    treatment = case_when(
      str_detect(sample_name, "AKO") ~ "cmAKO MI",
      str_detect(sample_name, "Control") ~ "WT Sham",
      str_detect(sample_name, "WT MI") ~ "WT MI"),
    region = case_when(
      str_detect(region, "border") ~ "Border Zone",
      str_detect(region, "infarct") ~ "Infarct Region",
      .default = region))

# Plot actinin area of all regions
remote.data <- data |> 
  select(sample_name, treatment, region, gene, area) |> 
  filter(gene == "Actinin")

# reorder norm type factor for plotting
remote.data <- remote.data |> 
  mutate(treatment = factor(treatment, levels = c("WT Sham", "WT MI", "cmAKO MI")),
         region = case_when(region == "remote" ~ "Remote Region",
                            .default = region))

# Plot for each gene
p.rem <- remote.data |> 
  ggplot(aes(x = treatment, y = area, fill = treatment, group = treatment)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", width = 0.9, color = "black") +
  geom_errorbar(aes(width = 0.9), stat='summary',  position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.7)  +
  facet_wrap(~region, ncol = 3) +
  labs(y = "Actinin area in remote region") +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#FF7F00")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        text = element_text(size = 11, color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(vjust = 0.5),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 9),
        strip.text = element_text(size = 11),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))

if(!dir.exists("results/supp_figs")){
  dir.create("results/supp_figs")
}
# Save plot to results 
png(file = "results/supp_figs/5a_actinin_rem.png",
    width = 8, 
    height = 3,
    units = "in",
    res = 800)

p.rem

dev.off()
