library(ggplot2)
library(PCAtools)
library(limma)

setwd(<working_directory>)

# Explore microarray dataset using Principal Components Analysis

# ------------------------ LOAD DATA ------------------------------------- #

# Load GSE15573 data set using GEOquery, here I loaded saved 
# data object from my system
load("datasets/GSE15573/GSE15573_ESet.RData")
show(eset)

# Select clinico-pathological variables from sample metadata
sample_data <- pData(eset)
sample_data <- data.frame( Trait = sample_data$characteristics_ch1,
                           Gender = sample_data$characteristics_ch1.1,
                           Age = sample_data$characteristics_ch1.2)

# Modify the variable definitions to make them more palatable
sample_data$Trait[sample_data$Trait == "status: Rheumatoid Arthritis Patient"] <- "RA"
sample_data$Trait[sample_data$Trait == "status: Control"] <- "Control"
sample_data$Gender <- gsub("gender: " , "", sample_data$Gender)
sample_data$Age <- gsub("age: " , "", sample_data$Age)
sample_data$Age <- as.numeric(sample_data$Age)

# Convert age to categorical variable
age_group<-cut(sample_data$Age, 
               breaks=c(44, 58, 81), right = T)
age_group
sample_data$age_group <- age_group

# Add modified variables to pheno data object
pData(eset)$Trait <- sample_data$Trait
pData(eset)$Gender <- sample_data$Gender
pData(eset)$Age <- sample_data$Age
pData(eset)$Age_group <- sample_data$age_group
head(pData(eset))

# -------------------------- PCA -------------------------------------------- #
# Run principal components analysis
expr <- exprs(eset)
expr <- normalizeVSN(expr) # Variance stabilization and log transform
metadata <- subset(pData(eset), select = c("Trait", "Gender", "Age_group")) 
p <- pca(expr, metadata = metadata, removeVar = 0.9)

# Create biplot
biplot(p, showLoadings = TRUE, lab = NULL)

# Create biplot and draw stat ellipses at 95% CI around groups
biplot(p,
       colby = 'Trait', colkey = c('Control' = 'forestgreen', 'RA' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)
ggsave("GSE15573_PCA_biplot.pdf", device = "pdf", units = "in", 
          width = 6, height = 6)

# Create pairs plot
pairsplot(p)
ggsave("GSE15573_pairs_plot.pdf", device = "pdf", units = "in", 
       width = 6, height = 6)


# Create loadings plot
plotloadings(p, labSize = 3)
ggsave("GSE15573_loadings_plot.pdf", device = "pdf", units = "in", 
       width = 6, height = 7)

# Plot loadings, fancy
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)
ggsave("GSE15573_loadings_plot_fancy.pdf", device = "pdf", units = "in", 
       width = 6, height = 7)

# Create eigencorrelation plot
pdf("GSE15573_eigencor_plot.pdf", width = 7, height = 7) 
eigencorplot(p,
             metavars = c('Trait','Gender','Age_group'))
dev.off()

# Identify optimal number of PCs
horn <- parallelPCA(expr) # Horn's method
horn$n # 7

elbow <- findElbowPoint(p$variance)
elbow # 7

# Create screeplot
screeplot(p, components = getComponents(p, 1:10),
          vline = c(horn$n, elbow)) + 
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))
ggsave("GSE15573_screeplot.pdf", device = "pdf", units = "in", 
       width = 7, height = 7)






