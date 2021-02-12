#################################################
#title: "TCGA_differential_methylation_analysis"
#author: "Nisha Rajasundaram"
#date: "Feb 10th 2021"
#output: CSV
################################################

# Clear workspace
rm(list=ls())

# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

#install.packages("DT")

# Load the required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)

# Search for all patients with DNA methylation (platform HumanMethylation450k) 
# and gene expression data (normal) for Colon Adenocarcinoma tumor (TCGA-COAD)

query.met <- GDCquery(project = "TCGA-COAD",
                      data.category = "DNA methylation",
                      data.type = "Methylation beta value",
                      sample.type = c("Primary Tumor"),
                      legacy = TRUE)

query.normal <- GDCquery(project = "TCGA-COAD",
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification", 
                      sample.type = "Solid Tissue Normal",
                      file.type  = "normalized_results",
                      legacy = TRUE)


# Get all patients that have DNA methylation and gene expression.
common.patients <- intersect(substr(getResults(query.met, cols = "cases"), 1, 12),
                             substr(getResults(query.normal, cols = "cases"), 1, 12))

# Only select the first 5 patients
query.met <- GDCquery(project = "TCGA-COAD",
                      data.category = "DNA methylation",
                      data.type = "Methylation beta value",
                      sample.type = c("Primary Tumor"),
                      legacy = TRUE,
                      barcode = common.patients[1:5])

query.normal <- GDCquery(project = "TCGA-COAD",
                         data.category = "Gene expression",
                         data.type = "Gene expression quantification", 
                         sample.type = "Solid Tissue Normal",
                         file.type  = "normalized_results",
                         legacy = TRUE,
                         barcode = common.patients[1:5])

# Checking for data
datatable(getResults(query.met, cols = c("data_type","cases")),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

datatable(getResults(query.normal, cols = c("data_type","cases")), 
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# Download queried files from TCGA
GDCdownload(query.met)
GDCdownload(query.normal)

# Pre-process the downloaded data
met.COAD <-GDCprepare(query.met)
normal.COAD <-GDCprepare(query.normal)

rowRanges(normal.COAD)
rowRanges(met.COAD)
colData(met.COAD)
colData(normal.COAD)

# Compartmentalize the downloaded data for further analysis 
data <- SummarizedExperiment(assays=SimpleList(counts=counts),
                             rowRanges=rowRanges, colData=colData)

# Differentially methylated regions Analysis using TCGAanalyze_DMC function
TCGAanalyze_DMC(
  data = data,
  groupCol = colData(data)
)

