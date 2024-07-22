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


# Plot for each gene
p.adj <- areas.adj |> 
  ggplot(aes(x = treatment, y = remote_normed, fill = norm_type, group = norm_type)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge", width = 0.9, color = "black") +
  geom_errorbar(aes(width = 0.9), stat='summary',  position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.7)  +
  facet_grid(gene ~ region  , scales = "free") +
  labs(y = "Ratio of spots compared to remote region") +
  #scale_fill_manual(values = c("#abd9e9","#fee090"), labels = c("Normalized to Actinin", "Actinin not considered")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        text = element_text(size = 8),
        legend.title = element_blank(),
        legend.position = c(0.25, 0.8),
        axis.title.y = element_text(size = 8))

if(!dir.exists("results/13_quant_images")){
  dir.create("results/13_quant_images")
}
# Save plot to results 
png(file = "results/13_quant_images/zbtb_norm.png",
    width = 4, 
    height = 5,
    units = "in",
    res = 800)

p.adj

dev.off()

areas |> 
  filter(gene == "Zbtb16") |> 
  ggplot(aes(x = treatment, y = value, fill = method)) +
  geom_col(position = "dodge") +
  facet_grid(region ~ gene, scales = "free") +
  labs(y = "Ratios",
       title = "Comparison of Pik3r1 and Zbtb16 expression by cardiac region during MI of cmAKO and WT mice",
       subtitle = "Mean number of RNAscope spots detected within the border and infarct zones were normalized to (divided by) the mean\n # of spots within the remote region within each tissue slice. This was repeated with a step to normalize the number of\nspots to the area of Actinin expression")


# For each ROI: 
# # of spots/actinin area / remote(# of spots/actinin area)
# # of spots/ remote(# of spots)