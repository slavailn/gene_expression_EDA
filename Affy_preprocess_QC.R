library(affy)
library(GEOquery)
library(ggplot2)
library(reshape2)
library(arrayQualityMetrics)
library(MatrixGenerics)

# Identify biomarkers to predict response to therapy 
# in polyarticular juvenile idiopathic arthritis (JIA) 
# using gene expression microarrays.
# 42 samples from 13 controls, 14 active patients, 9 patients in clinical 
# remission with medication (CRM), and 6 patients in clinical 
# remission without medication (CR). All patients had polyarticular JIA.


setwd(<working directory>)
setwd(<raw files>) # change to the directory containing .CEL files

# Download raw CEL files
getGEOSuppFiles("GSE15645")
# Untar and unzip CEL files from the command line

# Figure out which samples belong 2 which platform
gse <- getGEO("GSE15645")
pData(gse[[1]])$supplementary_file
write.csv(pData(gse[[1]]), file = "sample_info.csv")

# Add Trait and Disease_state variables
pData(gse[[1]])$Trait <- pData(gse[[1]])$`disease state:ch1`
pData(gse[[1]])$Disease_state <- pData(gse[[1]])$`symptom:ch1`

setwd("raw_files/")
AffyBatch <- ReadAffy()

# Run background correction, normalization and summarization
Affy_norm <- expresso(AffyBatch, bgcorrect.method="rma",
                        normalize.method="quantiles", pmcorrect.method="pmonly",
                        summary.method="medianpolish")
pData(Affy_norm) <- pData(gse[[1]])
# Save normalized ExpressionSet object 
save(Affy_norm, file="GSE15645_affy_norm.RData")

# Examine array qualities
# Create violin plot after normalization
# Get expression matrix
expr <- exprs(Affy_norm)
head(expr)
# Melt the expression matrix
expr_melt <- melt(as.data.frame(expr))
names(expr_melt) <- c("sample", "intensity")
head(expr_melt)

# The samples were background corrected but not normalized
ggplot(expr_melt, aes(x=sample, y = log2(intensity))) +
  geom_violin() + 
  geom_boxplot(width = 0.1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("GSE15645_intensities_violin_normalized.pdf", device = "pdf", units = "in", 
       width = 8, height = 5)

# Run array quality control analysis using arrayQualityMetrics
arrayQualityMetrics(Affy_norm,
                    outdir = "GSE15645_arrayQualityReport",
                    force = T,
                    do.logtransform = T,
                    intgroup = c("Disease_state"),
                    spatial = TRUE,
                    reporttitle = paste("arrayQualityMetrics report for", 
                                        deparse(substitute(expressionset)))
)

# ---------------------------------------------------------------------------- #
# Remove outlier and re-run the analysis
# Sample 14 is an outlier and must be removed
rownames(pData(gse[[1]]))[14]
# the offending sample - GSM391601

setwd("raw_files/")
AffyBatch <- ReadAffy()

# Run background correction, normalization and summarization
Affy_norm <- expresso(AffyBatch, bgcorrect.method="rma",
                      normalize.method="quantiles", pmcorrect.method="pmonly",
                      summary.method="medianpolish")
pheno <- pData(gse[[1]])
pheno <- pheno[-14,]
pData(Affy_norm) <- pheno
# Save normalized ExpressionSet object 
save(Affy_norm, file="GSE15645_affy_norm_outlier_removed.RData")

# Examine array qualities
# Create violin plot after normalization
# Get expression matrix
expr <- exprs(Affy_norm)
head(expr)

# Melt the expression matrix
expr_melt <- melt(as.data.frame(expr))
names(expr_melt) <- c("sample", "intensity")
head(expr_melt)

# The samples were background corrected but not normalized
ggplot(expr_melt, aes(x=sample, y = log2(intensity))) +
  geom_violin() + 
  geom_boxplot(width = 0.1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("GSE15645_intensities_violin_normalized_outlier_removed.pdf", 
       device = "pdf", units = "in", 
       width = 8, height = 5)

# Run array quality control analysis using arrayQualityMetrics
arrayQualityMetrics(Affy_norm,
                    outdir = "GSE15645_arrayQualityReport_outlier_removed",
                    force = T,
                    do.logtransform = T,
                    intgroup = c("Disease_state"),
                    spatial = TRUE,
                    reporttitle = paste("arrayQualityMetrics report for", 
                                        deparse(substitute(expressionset)))
)


