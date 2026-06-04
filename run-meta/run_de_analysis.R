library(dplyr)
library(limma)
source("run-meta/functions/meta_analysis_functions.R")


# Loading expression matrices
message("Loading expression data...")
dna_matrices <- readRDS("meta/matrices/dna_matrices.rds")
rna_matrices <- readRDS("meta/matrices/rna_matrices.rds") # this will contain expr and pData; use raw counts 'expr' for DESeq2 
message("Expression data loaded!")

# Loading pData
message("Loading pData...")
dna_pData <- readRDS("meta/pdata/dna_pData_with_condition.rds")
rna_pData <- readRDS("meta/pdata/rna_pData_with_condition.rds")
message("pData loaded!")

# Limma
message("Running DE analysis with limma...")
DE_res_dna <- run_limma_DE_list()
message("Analysis finished!")
