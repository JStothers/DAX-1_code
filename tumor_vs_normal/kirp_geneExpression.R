library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:/Users/jstot/Desktop/Jess_work/Bioinformatics/Project")

# queryKIRPpair<-GDCquery(project="TCGA-KIRP",
#                         legacy = TRUE, data.category = "Gene expression",
#                         data.type = "Gene expression quantification",
#                         platform = "Illumina HiSeq",
#                         sample.type = c("Primary Tumor",
#                                         "Solid Tissue Normal"),
#                         file.type = "results")
# 
# GDCdownload(queryKIRPpair)
# 
# prepKIRP <- GDCprepare(query = queryKIRPpair,
#                      save = TRUE, save.filename = "KIRP.rda",
#                      summarizedExperiment = TRUE)

load("KIRP.rda")
prepKIRP <- data

listall<-c(colnames(prepKIRP))
listall

#normal tissue, tumor primary
paired<-TCGAquery_MatchedCoupledSampleTypes(listall,c("NT","TP"))
paired

# querypaired<-GDCquery(project = "TCGA-KIRP",
#                       data.category = "Gene expression",
#                       data.type = "Gene expression quantification",
#                       experimental.strategy = "RNA-Seq",
#                       platform = "Illumina HiSeq",
#                       file.type = "results",
#                       barcode = paired,
#                       legacy = TRUE)
# 
# GDCdownload(queryKIRPpair)
# 
# KIRPpairedRNAseq<-GDCprepare(querypaired, save.filename = "KIRP_paired.rda")

load("KIRP_paired.rda")
KIRPpairedRNAseq <- data

KIRPpairedRNAseqmatrix<-assay(KIRPpairedRNAseq,"raw_count")
KIRPpairedRNAseqmatrix

dataNorm<-TCGAanalyze_Normalization(tabDF = KIRPpairedRNAseqmatrix,
                                    geneInfo = geneInfo)
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut = 0.25)
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
