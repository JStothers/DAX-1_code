# Case study from https://rdrr.io/bioc/TCGAbiolinks/f/vignettes/casestudy.Rmd

library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:/Users/jstot/Desktop/Jess_work/Bioinformatics/Project")
#-----------------------------------
# STEP 1: Search, download, prepare |
#-----------------------------------
# 1.1 - DNA methylation
# ----------------------------------
query.met <- GDCquery(project = "TCGA-ESCA", 
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450")
GDCdownload(query.met)


esca.met <- GDCprepare(query = query.met,
                      save = TRUE, 
                      save.filename = "escaDNAmet.rda",
                      summarizedExperiment = TRUE)

load("escaDNAmet.rda")
esca.met <- data
rm(data)

#-----------------------------------
# 1.2 - RNA expression
# ----------------------------------
query.exp <- GDCquery(project = "TCGA-ESCA", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results")
GDCdownload(query.exp)
esca.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "escaExp.rda")

load("escaExp.rda")
esca.exp <- data
rm(data)

# For DNA methylation, we perform a DMC (different methylated CpGs) analysis, which will give the difference 
# of DNA methylation for the probes of the groups and their significance value. The output can be seen in a volcano plot. 
# Note: Depending on the number of samples this function can be very slow due to the wilcoxon test, taking from hours to days.

# na.omit
esca.met <- subset(esca.met,subset = (rowSums(is.na(assay(esca.met))) == 0))

# Volcano plot
esca.met <- TCGAanalyze_DMC(esca.met, groupCol = "paper_DNA Methylation Cluster - EC",
                           group1 = "C1",
                           group2="C2",
                           p.cut = 10^-5,
                           diffmean.cut = 0.25,
                           legend = "State",
                           plot.filename = "MethylationCluster-EC_C1vsC2_metvolcano.png")

esca.met <- TCGAanalyze_DMC(esca.met, groupCol = "paper_SCNA High/Low",
                            group1 = "SCNA-High",
                            group2="SCNA-Low",
                            p.cut = 10^-5,
                            diffmean.cut = 0.25,
                            legend = "State",
                            plot.filename = "SCNA-HighvsSCNA-Low_metvolcano.png")
