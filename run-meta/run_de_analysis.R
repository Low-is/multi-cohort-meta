library(dplyr)
library(limma)
source("run-meta/functions/meta_analysis_functions.R")


# Loading expression matrices
dna_matrices <- readRDS("meta/matrices/dna_matrices.rds")
rna_matrices <- readRDS("meta/matrices/rna_matrices.rds") # this will contain expr and pData; use raw counts 'expr' for DESeq2 

# Loading pData


# Limma
DE_res_dna <- run_limma_DE_list()
