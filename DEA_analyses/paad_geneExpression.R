library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:/Users/jstot/Desktop/Jess_work/Bioinformatics/Project")

# queryPAADpair<-GDCquery(project="TCGA-PAAD",
#                         legacy = TRUE, data.category = "Gene expression",
#                         data.type = "Gene expression quantification",
#                         platform = "Illumina HiSeq",
#                         sample.type = c("Primary Tumor",
#                                         "Solid Tissue Normal"),
#                         file.type = "results")
# 
# GDCdownload(queryPAADpair)
# 
# prepPAAD <- GDCprepare(query = queryPAADpair,
#                      save = TRUE, save.filename = "PAAD.rda",
#                      summarizedExperiment = TRUE)

load("PAAD.rda")
prepPAAD <- data

listall<-c(colnames(prepPAAD))
listall

#normal tissue, tumor primary
paired<-TCGAquery_MatchedCoupledSampleTypes(listall,c("NT","TP"))
paired

# querypaired<-GDCquery(project = "TCGA-PAAD",
#                       data.category = "Gene expression",
#                       data.type = "Gene expression quantification",
#                       experimental.strategy = "RNA-Seq",
#                       platform = "Illumina HiSeq",
#                       file.type = "results",
#                       barcode = paired,
#                       legacy = TRUE)
# 
# GDCdownload(queryPAADpair)
# 
# PAADpairedRNAseq<-GDCprepare(querypaired, save.filename = "PAAD_paired.rda")

load("PAAD_paired.rda")
PAADpairedRNAseq <- data

PAADpairedRNAseqmatrix<-assay(PAADpairedRNAseq,"raw_count")
PAADpairedRNAseqmatrix

dataNorm<-TCGAanalyze_Normalization(tabDF = PAADpairedRNAseqmatrix,
                                    geneInfo = geneInfo)
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut = 0.20)
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("TP"))

dataDEGs<- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                           mat2 = dataFilt[,samplesTP],
                           Cond1type = "Normal",
                           Cond2type = "Tumor",
                           method = "glmLRT")

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],
                                          dataFilt[,samplesNT])

dataDEGsFiltLevel
dataDEGsFiltLevel["NR0B1",]
