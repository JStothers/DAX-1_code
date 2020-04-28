library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:/Users/jstot/Desktop/Jess_work/Bioinformatics/Project")

# queryLIHCpair<-GDCquery(project="TCGA-LIHC",
#                         legacy = TRUE, data.category = "Gene expression",
#                         data.type = "Gene expression quantification",
#                         platform = "Illumina HiSeq",
#                         sample.type = c("Primary Tumor",
#                                         "Solid Tissue Normal"),
#                         file.type = "results")
# 
# GDCdownload(queryLIHCpair)
# 
# prepLIHC <- GDCprepare(query = queryLIHCpair,
#                      save = TRUE, save.filename = "LIHC.rda",
#                      summarizedExperiment = TRUE)

load("LIHC.rda")
prepLIHC <- data

listall<-c(colnames(prepLIHC))
listall

#normal tissue, tumor primary
paired<-TCGAquery_MatchedCoupledSampleTypes(listall,c("NT","TP"))
paired

# querypaired<-GDCquery(project = "TCGA-LIHC",
#                       data.category = "Gene expression",
#                       data.type = "Gene expression quantification",
#                       experimental.strategy = "RNA-Seq",
#                       platform = "Illumina HiSeq",
#                       file.type = "results",
#                       barcode = paired,
#                       legacy = TRUE)
# 
# GDCdownload(queryLIHCpair)
# 
# LIHCpairedRNAseq<-GDCprepare(querypaired, save.filename = "LIHC_paired.rda")

load("LIHC_paired.rda")
LIHCpairedRNAseq <- data

LIHCpairedRNAseqmatrix<-assay(LIHCpairedRNAseq,"raw_count")
LIHCpairedRNAseqmatrix

dataNorm<-TCGAanalyze_Normalization(tabDF = LIHCpairedRNAseqmatrix,
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
