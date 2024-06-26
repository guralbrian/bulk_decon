# Load libs
library(remotes)
#install_version("ggplot2", version = "3.5.0")
libs <- c("ggtern", "tidyverse", "RColorBrewer", "reshape2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#### Loading and formatting of data ####
# Load compositions
decon.whole <- read.csv("data/processed/compositions/whole_samples.csv")
decon.frac <- read.csv("data/processed/compositions/fraction_samples.csv")

# refactorize the cell type levels
decon.frac <- decon.frac |> 
  mutate(cell.type = factor(cell.type))


#### Plot fractions ####
p.frac <- decon.frac |>
  ggplot(aes(x=new.id, y=Prop, fill=CellType))  +
  geom_bar(stat='identity',
           position = "fill",
           width = 1,
           color = "black")+
  scale_fill_brewer(name = "Cell Type",
                    palette = "Set2")  +
  facet_wrap(~cell.type,
             scales = "free_x")+
  ylab("Estimated Proportion") +
  theme_minimal() +
  theme(strip.text = element_text(size = 30),
        title = element_text(size = 20),
        axis.text.y  = element_text(size = 20),
        legend.text = element_text(size = 22),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right',
        legend.justification = 'center',
        legend.margin = margin(t = -5, r = 0, b = 0, l = 0, unit = "pt"),
        axis.text.x = element_blank(),
        plot.margin = margin(0, 5, 1, 5),
        plot.caption = element_text(size = 22, hjust = -0.3, face = "bold"),
        plot.tag = element_text(size = 22, hjust = -0.3),
        axis.title.y = element_text(size = 25),
        axis.title.x =  element_text(size = 25)) +
  xlab("Pure Cell Type Bulk RNAseq Replicates") 

png(file = "results/7_plot_comps/pure_cell_types.png",
    width = 16, 
    height = 8,
    units = "in",
    res = 600)

p.frac

dev.off()


#### Plot whole samples ####
#brewer.pal(n=8,"Paired")

my_palette <- c("#A6CEE3", "#1F78B4", "#FDBF6F", "#FF7F00")
legend.names <- c("Sham_1","Sham_2", "MI_1", "MI_2")


ordered_df <- decon.whole %>%
  filter(CellType == "Cardiomyocytes") %>%  # Filter rows where CellType is "Cardiomyocytes"
  arrange(Prop) %>%  # Sort these rows by Prop in descending order
  pull(new.id)  # Extract the new.id values in this order
# Use the ordered new.id values to reorder the original data frame
decon.whole <- decon.whole %>%
  mutate(new.id = factor(new.id, levels = ordered_df)) %>%
  arrange(new.id)
decon.whole$CellType_wrap = str_wrap(decon.whole$CellType, width = 12)


# Use the ordered new.id values to reorder the original data frame
decon.whole <- decon.whole %>%
  mutate(new.id = factor(new.id, levels = ordered_df)) %>%
  arrange(new.id) |> 
  mutate(Genotype_Treatment = paste(genotype, "-", treatment)) |> 
  mutate(Genotype_Treatment = factor(Genotype_Treatment, levels = c("WT - Sham", "WT - MI", "cmAKO - Sham", "cmAKO - MI")))
                                     
cell.type.order <-  decon.whole %>%
  group_by(CellType_wrap) |> 
  mutate(mean = mean(Prop)) |> 
  arrange(desc(mean)) |>
  pull(CellType_wrap) |> 
  unique()

design <- "AA
           BB"

# Generate boxplot
comp_celltype <- decon.whole   %>%
  ggplot(aes(x = factor(CellType_wrap, levels = as.character(cell.type.order)), y = Prop, fill = Genotype_Treatment)) +
  #geom_boxplot(position = position_dodge(0.9), width = 0.9, color = "black") +
  geom_bar(stat = "summary",position = position_dodge(0.9),  fun = mean,width = 0.9,  color = "black", alpha = 0.8) +
  geom_jitter(inherit.aes = T, 
              position = position_jitterdodge(jitter.width = 0.005, dodge.width = 0.9),
              size = 4, alpha = 0.5) +
  theme(axis.text.x = element_text(color = "black", angle =25, vjust = 0.5),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_blank(),
        legend.position = c(0.9, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), units = "cm"),
        text = element_text(size = 40)) +
  labs(y = "Estimated Proportion", 
       fill = "Treatment") +
  scale_fill_manual(values = my_palette)

# Find average comp at baseline
decon.whole |> 
  filter(Genotype_Treatment == "WT - Sham") |> 
  group_by(CellType) |> 
  summarize(mean = mean(Prop))

# Save plot to results 
#Save outputs
if(!dir.exists("results/7_plot_comps")){
  dir.create("results/7_plot_comps")
}

png(file = "results/7_plot_comps/sample_comps.png",
    width = 16, 
    height = 8,
    units = "in",
    res = 600)

comp_celltype

dev.off()


# format for ternary 
decon.wide <- decon.whole  |> 
  dplyr::select(new.id, CellType, Prop, Genotype_Treatment) |> 
  pivot_wider(names_from = "CellType", values_from = "Prop")


p.tern <- ggtern(decon.wide,aes(x=Cardiomyocytes,y=Fibroblast, z=Macrophage, color = Genotype_Treatment )) +
  geom_point(size = 10, alpha = 0.8) +
  scale_color_manual(values = my_palette) +
  theme_bw() +
  theme_hidegrid_minor() +
  theme_showarrows() +
  theme_hidetitles() +
  theme_arrowlarge() +
  tern_limit(T = 0.65,L = 1, R = 0.65) +
  theme(
    text = element_text(size = 28),
    legend.title = element_blank(),
    legend.position = "none",
    tern.axis.arrow = element_line(size = 3)
  )
png(file = "results/7_plot_comps/tern_whole.png",
    width = 8, 
    height = 8,
    units = "in",
    res = 600)

p.tern

dev.off()



# format for ternary 
decon.wide <- decon.frac  |> 
  dplyr::select(new.id, CellType, Prop, cell.type) |> 
  pivot_wider(names_from = "CellType", values_from = "Prop")

pal <- brewer.pal(3, "Set2")

p.tern.frac <- ggtern(decon.wide,aes(x=Cardiomyocytes,y=Fibroblast, z=`Endothelial Cells`, color = cell.type)) +
  geom_point(size = 13, alpha = 0.8) +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme_hidegrid_minor() +
  theme_showarrows() +
  theme_hidetitles() +
  theme_arrowlarge() +
  #tern_limit(T = 0.65,L = 1, R = 0.65) +
  theme(
    text = element_text(size = 28),
    legend.title = element_blank(),
    legend.position = "bottom",
    tern.axis.arrow = element_line(size = 3)
  )
png(file = "results/7_plot_comps/tern_frac.png",
    width = 8, 
    height = 8,
    units = "in",
    res = 600)

p.tern.frac

dev.off()

