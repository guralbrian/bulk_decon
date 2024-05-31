# Run dirichlet regression model and visualize the results

# Load libs and data 
libs <- c("tidyverse", "DirichletReg", "reshape2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)

#### Loading and formatting of data ####
# Load compositions
decon.whole <- read.csv("data/processed/compositions/whole_samples.csv") 



dir.model <- decon.whole |> 
  mutate(CellType = factor(CellType),
         Genotype = factor(genotype),
         Treatment = factor(treatment)
  ) 

## Prep DirichletReg matrix

# Add small value to remove zeros
#dir.model$Prop <- dir.model$Prop + 0.0001

# Make model matrix
dir.mat <- dir.model |>
  subset(select = c("new.id", "CellType", "Prop", "Treatment", "Genotype")) |>
  dcast(new.id + Treatment + Genotype ~ CellType, value.var = "Prop")

# Convert data to DirichletRegData object
dir.mat$CellTypes <- DR_data(dir.mat[,c(4:length(dir.mat))])

dir.mat$Treatment <- as.factor(dir.mat$Treatment) |>
  relevel(ref = "Sham") 

dir.mat$Genotype <- as.factor(dir.mat$Genotype) |>
  relevel(ref = "WT")

# Run Dirichlet regression
base <- 2
model.1 <- DirichReg(CellTypes ~ Treatment * Genotype, data = dir.mat, model = "alternative", base = base)
#model.2 <- DirichReg(CellTypes ~ Treatment * Genotype, data = dir.mat, model = "common")

# Save model coefficients 
summary(model.1) 
model.2 <- model.1
dir.results <- summary(model.2)[["coef.mat"]] |> as.data.frame(check.names = F)
rownames(dir.results) <- seq(1, nrow(dir.results))

# Clean the data frame
dir.results$Variable <- summary(model.2)[["coef.mat"]] |> rownames()
dir.results$CellType <- c(rep(summary(model.2)[["varnames"]][-base], each = 4),summary(model.2)[["varnames"]][base])
dir.results$Feature <- row.names(summary(model.2)[["coef.mat"]])
colnames(dir.results)[2] <- "StdError"

# Save
write.csv(dir.results, "data/processed/models/dirichelet_coefficients.csv", row.names = F)

# Prep labels for plot
dir.results <- dir.results |> 
  mutate(Feature_wrap = case_when(
    str_detect(Feature, "TreatmentMI$") ~ "MI",
    str_detect(Feature, ":") ~ "MI and\ncmAKO",
    str_detect(Feature, "GenotypecmAKO$") ~ "cmAKO",
    str_detect(Feature, "Intercept") ~ "Intercept"
  ))

# Plot non-intercept coefficients
#! Need to add significance labels
dir.results <- dir.results |> 
  subset(Feature != "(Intercept)")
err.plot <- dir.results |> 
  #subset(Feature != "(Intercept)" & `Pr(>|z|)` < 0.1) |> 
  ggplot(aes(y = Feature_wrap, x = Estimate, color = `Pr(>|z|)`, shape = CellType)) +
  geom_point(position = position_dodge(width = 0.6), size = 6) +
  geom_errorbarh(aes(xmin = Estimate - StdError, xmax = Estimate + StdError),
                 height = 0.2, position = position_dodge(width = 0.6), size = 2) +
  geom_errorbarh(data = dir.results[dir.results$unsig,], aes(xmin = Estimate - StdError, xmax = Estimate + StdError),
                 height = 0.2, position = position_dodge(width = 0.6), size = 2, color = "darkgrey") +
  scale_x_continuous(limits = c(-4, 2), n.breaks = 6) +
  theme_classic() +
  binned_scale(aesthetics = "color",
               scale_name = "stepsn", 
               palette = function(x) c("#FDE725FF", "#55C667FF", "#238A8DFF", "#482677FF", "grey"),
               breaks = c( 0.005, 0.01, 0.05, 0.1),
               limits = c(0.0005, 0.5),
               show.limits = T,
               guide = "colorsteps")+
  guides(shape = guide_legend(ncol = 1)) +
  #theme(legend.direction = "vertical", legend.box = "vertical")
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black", angle = 0),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        legend.position = "right",
        legend.justification = c("center", "center"),
        axis.title.y = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_line(colour = "grey", size = 0.5),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        text = element_text(size = 20)) +
  labs(x = "Coefficients (estimates of effect)", 
       shape = "Cell Types",
       color = "p-value") 

# Save 

#Save outputs
if(!dir.exists("results/8_dirichlet")){
  dir.create("results/8_dirichlet")
}

png(file = "results/8_dirichlet/dirichlet_coeff.png",
    width = 10, 
    height = 6,
    units = "in",
    res = 600)

err.plot

dev.off()



