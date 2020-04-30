# TCGAanalyze_DMR: Differentially methylated regions Analysis
# https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/analysis.html#tcgaanalyze_dmr:_differentially_methylated_regions_analysis

setwd("/Users/jstot/Desktop/Jess_work/Currently_Trying") 

library(TCGAbiolinks)

query <- GDCquery(project = "TCGA-GBM",
                  data.category = "DNA methylation",
                  platform = "Illumina Human Methylation 27",
                  legacy = TRUE, 
                  barcode = c("TCGA-02-0058-01A-01D-0186-05", "TCGA-12-1597-01B-01D-0915-05",
                              "TCGA-12-0829-01A-01D-0392-05", "TCGA-06-0155-01B-01D-0521-05",
                              "TCGA-02-0099-01A-01D-0199-05", "TCGA-19-4068-01A-01D-1228-05",
                              "TCGA-19-1788-01A-01D-0595-05", "TCGA-16-0848-01A-01D-0392-05"))
GDCdownload(query, method = "api")
data <- GDCprepare(query)

# "shortLetterCode" is a column in the SummarizedExperiment::colData(data) matrix
TCGAvisualize_meanMethylation(data, groupCol = "shortLetterCode",filename = NULL)

# setting y limits: lower 0, upper 1
TCGAvisualize_meanMethylation(data,groupCol = "shortLetterCode", 
                              filename = NULL, y.limits = c(0,1))

# setting y limits: lower 0
TCGAvisualize_meanMethylation(data,groupCol = "shortLetterCode", 
                              filename = NULL, y.limits = 0)

# Changing shapes of jitters to show subgroups
TCGAvisualize_meanMethylation(data,
                              groupCol = "paper_Pan.Glioma.DNA.Methylation.Cluster", 
                              subgroupCol ="vital_status", filename = NULL)


# Sorting bars by descending mean methylation
TCGAvisualize_meanMethylation(data,
                              groupCol = "paper_Pan.Glioma.DNA.Methylation.Cluster",
                              sort="mean.desc",
                              filename=NULL)

# Sorting bars by asc mean methylation
TCGAvisualize_meanMethylation(data,
                              groupCol = "paper_Pan.Glioma.DNA.Methylation.Cluster",
                              sort = "mean.asc",
                              filename=NULL)

TCGAvisualize_meanMethylation(data,
                              groupCol = "vital_status",
                              sort = "mean.asc",
                              filename=NULL)


