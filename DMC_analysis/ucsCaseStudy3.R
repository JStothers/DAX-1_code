# Case study from https://rdrr.io/bioc/TCGAbiolinks/f/vignettes/casestudy.Rmd

library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:/Users/jstot/Desktop/Jess_work/Bioinformatics/Project")
#-----------------------------------
# STEP 1: Search, download, prepare |
#-----------------------------------
# 1.1 - DNA methylation
# ----------------------------------
query.met <- GDCquery(project = "TCGA-UCS", 
                      legacy = TRUE,
                      data.category = "DNA methylation",
                      platform = "Illumina Human Methylation 450")
GDCdownload(query.met)

ucs.met <- GDCprepare(query = query.met,
                      save = TRUE, 
                      save.filename = "ucsDNAmet.rda",
                      summarizedExperiment = TRUE)

#-----------------------------------
# 1.2 - RNA expression
# ----------------------------------
query.exp <- GDCquery(project = "TCGA-UCS", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results")
GDCdownload(query.exp)
ucs.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "ucsExp.rda")

# For DNA methylation, we perform a DMC (different methylated CpGs) analysis, which will give the difference 
# of DNA methylation for the probes of the groups and their significance value. The output can be seen in a volcano plot. 
# Note: Depending on the number of samples this function can be very slow due to the wilcoxon test, taking from hours to days.

# na.omit
ucs.met <- subset(ucs.met,subset = (rowSums(is.na(assay(ucs.met))) == 0))

# Volcano plot
ucs.met <- TCGAanalyze_DMC(ucs.met, groupCol = "paper_MethylationCluster",
                           group1 = "1",
                           group2="2",
                           p.cut = 10^-5,
                           diffmean.cut = 0.25,
                           legend = "State",
                           plot.filename = "MethyCluster1vs2_metvolcano.png")

ucs.met <- TCGAanalyze_DMC(ucs.met, groupCol = "paper_MethylationCluster",
                           group1 = "1",
                           group2="3",
                           p.cut = 10^-5,
                           diffmean.cut = 0.25,
                           legend = "State",
                           plot.filename = "MethyCluster1vs3_metvolcano.png")

ucs.met <- TCGAanalyze_DMC(ucs.met, groupCol = "paper_MethylationCluster",
                           group1 = "2",
                           group2="3",
                           p.cut = 10^-5,
                           diffmean.cut = 0.25,
                           legend = "State",
                           plot.filename = "MethyCluster2vs3_metvolcano.png")
