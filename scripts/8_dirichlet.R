# Run dirichlet regression model and visualize the results

# Load libs and data 

# Load libs
libs <- c("tidyverse", "DirichletReg", "reshape2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#### Loading and formatting of data ####
# Load compositions
decon.whole <- read.csv("data/processed/compositions/whole_samples.csv")


dir.model <- decon.whole |> 
  mutate(CellType = factor(CellType),
         Genotype = factor(Genotype),
         Treatment = factor(Treatment)
  ) 

## Prep DirichletReg matrix

# Add small value to remove zeros
dir.model$Prop <- dir.model$Prop + 0.0001

# Make model matrix
dir.mat <- dir.model |>
  subset(select = c("Sub", "CellType", "Prop", "Treatment", "Genotype")) |>
  dcast(Sub + Treatment + Genotype ~ CellType, value.var = "Prop")


# Convert data to DirichletRegData object
dir.mat$CellTypes <- DR_data(dir.mat[,c(4:length(dir.mat))])

dir.mat$Treatment <- as.factor(dir.mat$Treatment) |>
  relevel(ref = "Sham") 

dir.mat$Genotype <- as.factor(dir.mat$Genotype) |>
  relevel(ref = "WT")
# Run Dirichlet regression
#model.1 <- DirichReg(CellTypes ~ Treatment * Genotype, data = dir.mat, model = "alternative", base = 3)
model.2 <- DirichReg(CellTypes ~ Treatment * Genotype, data = dir.mat, model = "common")

# Save model coefficients 
summary(model.2) # Model results don't look great, might want to remove outlier
dir.results <- summary(model.2)[["coef.mat"]] |> as.data.frame(check.names = F)
rownames(dir.results) <- seq(1, nrow(dir.results))



# Clean the data frame
dir.results$Variable <- summary(model.2)[["coef.mat"]] |> rownames()
dir.results$CellType <- rep(summary(model.2)[["varnames"]], each = 4)
dir.results$Feature <- row.names(summary(model.2)[["coef.mat"]])
colnames(dir.results)[2] <- "StdError"

# Save
write.csv(dir.results, "data/processed/models/dirichelet_coefficients.csv", row.names = F)

# Prep labels for plot
dir.results <- dir.results |> 
  mutate(Feature_wrap = case_when(
    str_detect(Feature, "TreatmentTAC$") ~ "Myocardial\nInfarction",
    str_detect(Feature, ":") ~ "Myocardial\nInfarction\nand KO",
    str_detect(Feature, "GenotypeKO$") ~ "A1-AR KO",
    str_detect(Feature, "Intercept") ~ "Intercept"
  ))



# Plot non-intercept coefficients
#! Need to add significance labels
err.plot <- dir.results |> 
  subset(Feature != "(Intercept)") |> #& `Pr(>|z|)` < 0.05) |> 
  ggplot(aes(y = Feature_wrap, x = Estimate, color = CellType)) +
  geom_point(position = position_dodge(width = 0.6), size = 6) +
  geom_errorbarh(aes(xmin = Estimate - StdError, xmax = Estimate + StdError),
                 height = 0.2, position = position_dodge(width = 0.6), size = 2) +
  # geom_text(aes(label = Significance, x = Estimate + StdError, group = CellType), 
  #          position = position_dodge(width = 0.6), 
  #          colour = "black", hjust = -0.2, size = 10) +
  scale_color_brewer(name = "Cell Type",
                     palette = "Dark2")  +
  scale_x_continuous(limits = c(-5, 3), n.breaks = 8) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20, angle = 0),
        axis.ticks = element_blank(),
        legend.position = c(0.05,0.2),
        legend.justification = c("left", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.title.y = element_text(color = "black", size = 24),
        axis.title.x = element_text(color = "black", size = 24),
        #panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(colour = "grey", size = 0.5),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='white')) +
  labs(x = "Coeffiecients (estimates of effect)", 
       y = "Predictor Variables", 
       legend = "Cell Types") 

# Save 
png(file = "results/8_dirichlet/dirichlet_coeff.png",
    width = 1200, 
    height = 800,
    units = "px",
    res = 100)

err.plot

dev.off()



