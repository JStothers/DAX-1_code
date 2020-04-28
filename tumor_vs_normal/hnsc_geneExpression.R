library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:/Users/jstot/Desktop/Jess_work/Bioinformatics/Project")

# queryHNSCpair<-GDCquery(project="TCGA-HNSC",
#                         legacy = TRUE, data.category = "Gene expression",
#                         data.type = "Gene expression quantification",
#                         platform = "Illumina HiSeq",
#                         sample.type = c("Primary Tumor",
#                                         "Solid Tissue Normal"),
#                         file.type = "results")
# 
# GDCdownload(queryHNSCpair)
# 
# prepHNSC <- GDCprepare(query = queryHNSCpair,
#                      save = TRUE, save.filename = "HNSC.rda",
#                      summarizedExperiment = TRUE)

load("HNSC.rda")
prepHNSC <- data

listall<-c(colnames(prepHNSC))
listall

#normal tissue, tumor primary
paired<-TCGAquery_MatchedCoupledSampleTypes(listall,c("NT","TP"))
paired

# querypaired<-GDCquery(project = "TCGA-HNSC",
#                       data.category = "Gene expression",
#                       data.type = "Gene expression quantification",
#                       experimental.strategy = "RNA-Seq",
#                       platform = "Illumina HiSeq",
#                       file.type = "results",
#                       barcode = paired,
#                       legacy = TRUE)
# 
# GDCdownload(queryHNSCpair)
# 
# HNSCpairedRNAseq<-GDCprepare(querypaired, save.filename = "HNSC_paired.rda")

load("HNSC_paired.rda")
HNSCpairedRNAseq <- data

HNSCpairedRNAseqmatrix<-assay(HNSCpairedRNAseq,"raw_count")
HNSCpairedRNAseqmatrix

dataNorm<-TCGAanalyze_Normalization(tabDF = HNSCpairedRNAseqmatrix,
                                    geneInfo = geneInfo)
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut = 0.15)
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

