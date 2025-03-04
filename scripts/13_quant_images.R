# This script is meant to process the table from the RNAscope image quantification
# For each sample, three zones were considered: border, remote, and infarct
# 4-5 reprentative areas were chosen within each zone and the number of distinct spots were tallied

# Load data
data <- read.csv("data/raw/rnascope/alpha1_ar_image_quant.csv")

# Load libraries
library(tidyverse)

# Pull gene, sample type, and number from sample name

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

# For each ROI,region,treatment, get gene spots normed to actinin area

cm.normed <- data |> 
  filter(gene.group == "Actinin") |> 
  select(treatment, region, ROI, gene, area, number) |> 
  pivot_wider(values_from = c("area", "number"), names_from = gene) |> 
  mutate(cm_Zbtb16 = number_Zbtb16/area_Actinin,
         cm_Pik3r1 = number_Pik3r1/area_Actinin) |> 
  select(-area_Actinin, -area_Pik3r1, -area_Zbtb16, -number_Actinin) |> 
  pivot_longer(cols = starts_with("number_"), names_prefix = "number_", names_to = "gene", values_to = "counts") |> 
  pivot_longer(cols = starts_with("cm_"), names_prefix = "cm_", names_to = "gene2", values_to = "cm_norm_counts") |> 
  filter(gene == gene2) |> # Fix weird pivoting
  select(-gene2)  |> 
  pivot_longer(cols = contains("counts"), names_to = "norm_type") |> 
  mutate(
    norm_type = case_when(
      norm_type == "counts" ~ "raw_counts",
      norm_type == "cm_norm_counts" ~ "cm_norm_counts"
    )
  )


areas.adj <- cm.normed |> 
  filter(region == "remote") |> 
  group_by(treatment, region, gene, norm_type) |> 
  summarize(remote_value = mean(value)) |>
  ungroup() |> 
  select(-region) |> 
  right_join(cm.normed, by = c("treatment", "gene", "norm_type")) |> 
  filter(region != "remote") |> 
  mutate(remote_normed = value / remote_value) |> 
  select(-remote_value, -value) 

# take ratio of CMs (channel 2 actinin) to number of dots detected for channel 3/4, for each ROI
# Then, take ratio of infarct/border zones to remote regions and compare across sample.types

# reorder norm type factor for plotting
areas.adj <- areas.adj |> 
  mutate(norm_type = factor(norm_type, levels = c("raw_counts", "cm_norm_counts")),
         treatment = factor(treatment, levels = c("WT Sham", "WT MI", "cmAKO MI")))

# Plot for each gene
p.adj <- areas.adj |> 
  filter(region == "Infarct Region") |> 
  ggplot(aes(x = treatment, y = remote_normed, fill = norm_type, group = norm_type)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", width = 0.9, color = "black") +
  geom_errorbar(aes(width = 0.9), stat='summary',  position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.7)  +
  facet_wrap( ~  gene , scales = "free") +
  labs(y = "INF:remote-region RNA spot counts") +
  scale_fill_manual(values = c("#fee090", "#abd9e9"), labels = c("Actinin not considered", "Normalized to Actinin")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        text = element_text(size = 11, color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(vjust = 0.5),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 9),
        strip.text = element_text(size = 11),
        legend.position = "bottom",
        legend.justification = c("center", "center"),
        legend.box.just = "center",
        legend.text = element_text(size = 12, margin = margin(t = 0, r = 5, b = 0, l = 5, unit = "pt")),
        legend.key.size = unit(0.8, "lines"),  # Adjust size of the keys
        legend.spacing.x = unit(1, "mm"),    # Reduce spacing between columns of the legend
        legend.spacing.y = unit(0.5, "mm"),    # Reduce spacing between rows of the legend
        legend.box.spacing = unit(0.5, "mm"),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))

if(!dir.exists("results/13_quant_images")){
  dir.create("results/13_quant_images")
}
# Save plot to results 
png(file = "results/13_quant_images/zbtb_norm_inf.png",
    width = 8, 
    height = 3,
    units = "in",
    res = 800)

p.adj

dev.off()


# Save the plot with a transparent background
ggsave("results/13_quant_images/gmb_2024_zbtb_norm_inf.png",
       width = 13.81/3,
       height = 7.37/3,
       units = "in",
       dpi = 600,
       plot = p.adj, bg = "transparent")

# Plot for each gene
p.adj <- areas.adj |> 
  filter(region == "Border Zone") |> 
  ggplot(aes(x = treatment, y = remote_normed, fill = norm_type, group = norm_type)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", width = 0.9, color = "black") +
  geom_errorbar(aes(width = 0.9), stat='summary',  position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.7)  +
  facet_wrap( ~  gene , scales = "free") +
  labs(y = "Ratio of BZ and remote region RNA spot counts") +
  scale_fill_manual(values = c("#fee090", "#abd9e9"), labels = c("Actinin not considered", "Normalized to Actinin")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        text = element_text(size = 11, color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(vjust = 0.5),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 9),
        strip.text = element_text(size = 11),
        legend.position = "bottom",
        legend.justification = c("center", "center"),
        legend.box.just = "center",
        legend.text = element_text(size = 12, margin = margin(t = 0, r = 5, b = 0, l = 5, unit = "pt")),
        legend.key.size = unit(0.8, "lines"),  # Adjust size of the keys
        legend.spacing.x = unit(1, "mm"),    # Reduce spacing between columns of the legend
        legend.spacing.y = unit(0.5, "mm"),    # Reduce spacing between rows of the legend
        legend.box.spacing = unit(0.5, "mm"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(color = "darkgrey"),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), units = "cm"))

if(!dir.exists("results/supp_figs")){
  dir.create("results/supp_figs")
}
# Save plot to results 
png(file = "results/supp_figs/zbtb_norm_bz.png",
    width = 8, 
    height = 3,
    units = "in",
    res = 800)

p.adj

dev.off()

