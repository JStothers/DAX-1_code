# based on Sue and Rashmy's code

#set the working directory
setwd("C:/Users/jstot/Desktop/Jess_work/Bioinformatics/Project")


#install.packages("devtools")
# devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks",dependencies = T) 
library(TCGAbiolinks)
library(DESeq2)
library(GenomicRanges)
library(DESeq)


#Gene expression aligned against hg19.

#query.exp.hg19 <- GDCquery(project = "TCGA-BRCA",
#                           data.category = "Gene expression",
#                           data.type = "Gene expression quantification",
#                           file.type  = "results",
#                           experimental.strategy = "RNA-Seq",
#                           sample.type = c("Primary solid Tumor","Solid Tissue Normal"),
#                            legacy = TRUE)
#GDCdownload(query.exp.hg19)
#TCGA_raw <- GDCprepare(query.exp.hg19, save = TRUE, save.filename = "TCGA_raw.rda") # yeilds a SummarizedExperiment object

load("TCGA_raw.rda") # the data has info about all the samples i.e. 
# "Primary solid Tumor","Solid Tissue Normal", "metastatic" 


# get assay matrix information
TCGA_raw <- data
TCGA_raw.matrix <- assay(TCGA_raw) 
barcodes_RNASeq_GDC <- colnames(TCGA_raw.matrix) 

#download clinical data
clinical_download <- function(barcodes){
  if(!is.null(barcodes) && length(barcodes)>0){
    if( nchar(barcodes[1])>12){
      barcodes <- substr(barcodes,1,12)
    }
    clin.query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", 
                           barcode = barcodes)
    tryCatch(GDCdownload(clin.query), error = function(e) 
      GDCdownload(clin.query, method = "client"))
    clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
    #er_pos <- TCGAquery_clinicFilt( barcodes, clinical.patient, ER="Positive")
    #er_neg <- TCGAquery_clinicFilt( barcodes, clinical.patient, ER="Negative")
    return (clinical.patient)
  }
}

# Downloads clinical data using the barcodes and filters ER positive and negative
clinical_data <- clinical_download(barcodes_RNASeq_GDC)

# get sample information
# sample_data <- as.data.frame(colData(TCGA_raw)) 
sample_data <- colData(TCGA_raw) 

# get feature matrix information (gene information)
gene_data <- rowRanges(TCGA_raw) 

# converting gene_data , a GRangesList object to a data frame
gdata <- as.data.frame(gene_data, row.names = NULL, optional = FALSE, 
                       value.name = "value", use.outer.mcols = FALSE, 
                       group_name.as.factor = FALSE)

# Use "sample_data$subtype_PAM50.mRNA" to extract data based on molecular subtypes.
# levels(sample_data$subtype_PAM50.mRNA)
# [1] "Basal-like"    "HER2-enriched" "Luminal A"     "Luminal B"     "Normal-like"

# Subset data to perform DEA and extract barcodes for different subtypes for further
# processing
basal_like_data <- subset(sample_data, paper_BRCA_Subtype_PAM50 == "Basal")
barcodes_basal <- rownames(basal_like_data)

non_basal <- subset(sample_data, paper_BRCA_Subtype_PAM50 != "Basal" | is.na(paper_BRCA_Subtype_PAM50))
barcodes_other <- rownames(non_basal)

HER2_data <- subset(sample_data, paper_BRCA_Subtype_PAM50 == "Her2")
barcodes_Her2 <- rownames(HER2_data)

lumi_A_data <- subset(sample_data, paper_BRCA_Subtype_PAM50 == "LumA")
barcodes_lum_A <- rownames(lumi_A_data)

lumi_B_data <- subset(sample_data, paper_BRCA_Subtype_PAM50 == "LumB")
barcodes_lum_B <- rownames(lumi_B_data)

norm_like_data <- subset(sample_data, paper_BRCA_Subtype_PAM50 == "Normal")
barcodes_normal <- rownames(norm_like_data)

barcodes_na <- rownames(subset(sample_data, is.na(paper_BRCA_Subtype_PAM50)))



TCGAvisualize_meanMethylation(data, groupCol = "paper_BRCA_Subtype_PAM50",filename = NULL, title = "Mean DNA Methylation", 
                              ylab = expression(paste("Mean DNA Methylation (", beta, "-values)")))



