
setwd("/Users/jstot/Desktop/Jess_work")


library(GenomicRanges)
library(DESeq)
library(DESeq2)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(stringi)
library(hms)
library(ELMER.data)
library(ELMER)
library(parallel)
library(DT)
library(dplyr)
library(BiocStyle)
library(GenomicInteractions)
library(rtracklayer)
library(AnnotationDbi)
library(GenomicFeatures)

#########search , download,prepare. Once we saved it as methy.rda we do not need to run this code anymore. We just load the data#########################
####DNA methylation#############
# mydir='/Users/jstot/Desktop/Jess_work'
# query <- GDCquery(project = "TCGA-BRCA",
#                    data.category = "DNA methylation",
#                    platform = "Illumina Human Methylation 27",
#                    legacy = TRUE)
# GDCdownload(query, directory = mydir)
# methy_data <- GDCprepare(query, save = TRUE,save.filename = "methy.rda",
#                           directory = mydir,
#                           remove.files.prepared = TRUE)

load("methy.rda")
Meth_data<- data

colMAE <- colData(Meth_data)
Gene_data <- rowRanges(data)

sample_gene <- data.frame(Gene_data)
sample_col <- as.data.frame(colMAE)
sample_elementmeta <- as.data.frame(elementMetadata(data))
sample_metadata <- as.data.frame(metadata(data))
sample_meth__data <- as.data.frame(assays(data))


TCGAvisualize_meanMethylation(data=Meth_data,
                              groupCol = "tumor_stage", #Columns that defines groups, must be from colData
                              title = "Mean DNA Methylation by tumor Stage",
                              filename =NULL,
                              ylab = expression(paste("Mean DNA Methylation (", beta, "-values)")))


TCGAvisualize_meanMethylation(data=Meth_data,
                              groupCol = "definition", #Columns that defines groups
                              title = "Mean DNA Methylation by Definition",
                              filename =NULL,
                              ylab = expression(paste("Mean DNA Methylation (", beta, "-values)")))



TCGAvisualize_meanMethylation(data=Meth_data,
                              groupCol = "shortLetterCode", #Columns that defines groups
                              title = "Mean DNA Methylation by Shortlettercode",
                              filename =NULL,
                              ylab = expression(paste("Mean DNA Methylation (", beta, "-values)")))







