# Case study from https://rdrr.io/bioc/TCGAbiolinks/f/vignettes/casestudy.Rmd

library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:/Users/jstot/Desktop/Jess_work/Bioinformatics/Project")
#-----------------------------------
# STEP 1: Search, download, prepare |
#-----------------------------------
# 1.1 - DNA methylation
# ----------------------------------
query.met <- GDCquery(project = "TCGA-TGCT", 
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450")
GDCdownload(query.met)


tgct.met <- GDCprepare(query = query.met,
                      save = TRUE, 
                      save.filename = "tgctDNAmet.rda",
                      summarizedExperiment = TRUE)

load("tgctDNAmet.rda")
tgct.met <- data
rm(data)

#-----------------------------------
# 1.2 - RNA expression
# ----------------------------------
query.exp <- GDCquery(project = "TCGA-TGCT", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results")
GDCdownload(query.exp)
tgct.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "tgctExp.rda")

load("tgctExp.rda")
tgct.exp <- data
rm(data)

# For DNA methylation, we perform a DMC (different methylated CpGs) analysis, which will give the difference 
# of DNA methylation for the probes of the groups and their significance value. The output can be seen in a volcano plot. 
# Note: Depending on the number of samples this function can be very slow due to the wilcoxon test, taking from hours to days.

# na.omit
tgct.met <- subset(tgct.met,subset = (rowSums(is.na(assay(tgct.met))) == 0))

# Volcano plot
tgct.met <- TCGAanalyze_DMC(tgct.met, groupCol = "paper_MethyLevel",
                           group1 = "CIMP-high",
                           group2="CIMP-low",
                           p.cut = 10^-5,
                           diffmean.cut = 0.25,
                           legend = "State",
                           plot.filename = "CIMP-highvsCIMP-low_metvolcano.png")

