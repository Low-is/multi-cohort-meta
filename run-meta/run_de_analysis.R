library(dplyr)
library(limma)
library(yaml)
library(AnnotationDbi)
library(org.Hs.eg.db)
source("run-meta/functions/meta_analysis_functions.R")


# Loading config file
config <- read_yaml("run-meta/config/config.yaml")
gene_config <- config$limma$genes

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

# Loading meta genes
meta_genes <- readRDS(gene_config$meta_genes)
# Or if leaving genes_of_interest NULL
all_genes <- gene_config$all

names(dna_pData)

# Limma
#message("Running DE analysis with limma...")
#DE_res_dna <- run_limma_DE_list(dna_matrices,
                               dna_pData,
                               genes_of_interest = all_genes)
#message("Analysis finished!")
#message("Saving limma results...")
#write.csv(DE_res_dna, "run-meta/output/limma_res.csv")
