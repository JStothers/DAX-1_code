library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:/Users/jstot/Desktop/Jess_work/Bioinformatics/Project")

# queryBLCApair<-GDCquery(project="TCGA-BLCA",
#                         legacy = TRUE, data.category = "Gene expression",
#                         data.type = "Gene expression quantification",
#                         platform = "Illumina HiSeq",
#                         sample.type = c("Primary Tumor",
#                                         "Solid Tissue Normal"),
#                         file.type = "results")
# 
# GDCdownload(queryBLCApair)
# 
# prepBLCA <- GDCprepare(query = queryBLCApair,
#                      save = TRUE, save.filename = "BLCA.rda",
#                      summarizedExperiment = TRUE)

load("BLCA.rda")
prepBLCA <- data

listall<-c(colnames(prepBLCA))
listall

#normal tissue, tumor primary
paired<-TCGAquery_MatchedCoupledSampleTypes(listall,c("NT","TP"))
paired

# querypaired<-GDCquery(project = "TCGA-BLCA",
#                       data.category = "Gene expression",
#                       data.type = "Gene expression quantification",
#                       experimental.strategy = "RNA-Seq",
#                       platform = "Illumina HiSeq",
#                       file.type = "results",
#                       barcode = paired,
#                       legacy = TRUE)
# 
# GDCdownload(queryBLCApair)
# 
# BLCApairedRNAseq<-GDCprepare(querypaired, save.filename = "BLCA_paired.rda")

load("BLCA_paired.rda")
BLCApairedRNAseq <- data

BLCApairedRNAseqmatrix<-assay(BLCApairedRNAseq,"raw_count")
BLCApairedRNAseqmatrix

dataNorm<-TCGAanalyze_Normalization(tabDF = BLCApairedRNAseqmatrix,
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
