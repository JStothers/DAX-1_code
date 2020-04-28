library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:/Users/jstot/Desktop/Jess_work/Bioinformatics/Project")
#-----------------------------------
# STEP 1: Search, download, prepare |
#-----------------------------------
# 1.1 - DNA methylation
# ----------------------------------
query.met <- GDCquery(project = "TCGA-SARC", 
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450")
GDCdownload(query.met)


sarc.met <- GDCprepare(query = query.met,
                      save = TRUE, 
                      save.filename = "sarcDNAmet.rda",
                      summarizedExperiment = TRUE)

load("sarcDNAmet.rda")
sarc.met <- data
rm(data)

#-----------------------------------
# 1.2 - RNA expression
# ----------------------------------
query.exp <- GDCquery(project = "TCGA-SARC", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results")
GDCdownload(query.exp)
sarc.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "sarcExp.rda")

# For DNA methylation, we perform a DMC (different methylated CpGs) analysis, which will give the difference 
# of DNA methylation for the probes of the groups and their significance value. The output can be seen in a volcano plot. 
# Note: Depending on the number of samples this function can be very slow due to the wilcoxon test, taking from hours to days.

# na.omit
sarc.met <- subset(sarc.met,subset = (rowSums(is.na(assay(sarc.met))) == 0))

# Volcano plot
sarc.met <- TCGAanalyze_DMC(sarc.met, groupCol = "paper_MethyLevel",
                           group1 = "CIMP-high",
                           group2="CIMP-low",
                           p.cut = 10^-5,
                           diffmean.cut = 0.25,
                           legend = "State",
                           plot.filename = "CIMP-highvsCIMP-low_metvolcano.png")
