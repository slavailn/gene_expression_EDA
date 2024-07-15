library(GEOquery)
setwd(<working_directory>)

# Download existing RA datasets
# GEO accession: GSE93776 
#   
# Status	Public on Jul 16, 2018
# Title	Immune cells gene expression from rheumatoid arthritis and healthy donors
# Experiment type	Expression profiling by array
# Summary	Sustained clinical remission (CR) without drug 
# treatment has not been achieved in patients with rheumatoid arthritis (RA). 
# This implies a substantial difference between CR and the healthy state, 
# but it has yet to be quantified. We report a longitudinal monitoring of the 
# drug response at multi-omics levels in the peripheral blood of patients with RA. 
# Our data reveal that drug treatments alter the molecular profile closer to that
# of HCs at the transcriptome, serum proteome and immunophenotype level. 
# Patient follow-up suggests that the molecular profile after drug treatments 
# is associated with long-term stable CR. 


# Cite: Tasaki S, Suzuki K, Kassai Y, Takeshita M et al. Multi-omics 
# monitoring of drug response in rheumatoid arthritis in pursuit of molecular 
# remission. Nat Commun 2018 Jul 16;9(1):2755. PMID: 30013029

# Download GSE93776
gse <- getGEO("GSE93776")
show(gse)
gse[[1]]
show(pData(phenoData(gse[[1]]))[1:4,])

# Get sample data
eset <- gse[[1]]
write.csv(pData(eset), file="GSE93776_phenodata.csv")
save(eset, file="GSE93776_ESet_no_norm.RData")

# Get annotation
gpl <- getGEO(eset@annotation)
class(gpl)
head(gpl@dataTable@table)
annotation <- gpl@dataTable@table
# Save annotation
write.csv(annotation, file="GSE93776_annotation.csv")

# Download raw files
getGEOSuppFiles("GSE93776", makeDirectory = TRUE, baseDir = getwd(),
                fetch_files = TRUE, filter_regex = NULL)



