library(yaml)
library(jsonlite)
library(COCONUT)
library(DESeq2)
library(data.table)
library(dplyr)
source("run-meta/functions/meta_analysis_functions.R")
source("run-meta/functions/Sepsis_MC_analysis_functions.R")
source("meta/functions/expr_mtx.R")
source("meta/functions/filter_data.R")


# Loading config file
message("Loading config file...")
config <- yaml::read_yaml("run-meta/config/config.yaml")
message("Config file loaded!")

# Loading named list of studies
message("Loading list of studies...")
dna_studies <- jsonlite::fromJSON(config$analysis$input$dna_gse_file)
rna_studies <- jsonlite::fromJSON(config$analysis$input$rna_gse_file)
message("Studies loaded!")

# Loading expression matrices
message("Loading expression matrices...")
dna_matrices <-  readRDS("meta/matrices/dna_matrices.rds")
dna_matrices <- dna_matrices[!sapply(dna_matrices, is.null)]

# Function to get rid of methylation experiements
is_meth <- function(genes) {
  mean(grepl("^cg", genes, ignore.case = TRUE)) > 0.05
}

dna_matrices <- dna_matrices[!sapply(dna_matrices, function(m) is_meth(rownames(m)))]

#rna_matrices <- readRDS("meta/matrices/rna_matrices.rds") # the str is a list of 2: x$expr and x$pData
#rna_matrices <- rna_matrices[!sapply(rna_matrices, is.null)] # Removing NULL single cell datasets
rna_matrices <- readRDS("meta/matrices/norm_rna_mtxs.rds")
message("Matrices loaded!")

                                     
message("Loading pData...")
dna_pData <- readRDS("meta/pdata/dna_pData_with_condition.rds")
rna_pData <- readRDS("meta/pdata/rna_pData_with_condition.rds")
message("pData loaded!")


# Need to add code that filters matrices to match dimmensions of pData         
# Filtering list
dna_matrices <- dna_matrices[names(dna_pData)] 
rna_matrices <- rna_matrices[names(rna_pData)]


# ----------------------------
# NORMALIZATION (GLOBAL SAFE)
# ----------------------------
normalize <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9 ]", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

#generate_candidate_titles(rna_pData[["GSE106910"]]$title)
colnames(rna_matrices[["GSE106910"]])
normalize(colnames(rna_matrices[["GSE106910"]]))


                                                                                                   
#message("Getting norm RNA counts...")
#norm_rna_mtxs <- get_norm_RNA_counts(rna_matrices, pData = rna_pData)
#message("Extracted norm RNA counts!")


# Find common genes across all studies being used for meta-analysis
#message("Searching for common genes...")
#common_genes <- find_common_genes(DNA = config$analysis$modalities$DNA,
                                  #RNA = config$analysis$modalities$RNA,
                                  #list_of_dna_mtx = dna_matrices,
                                  #list_of_rna_mtx = norm_rna_mtxs,
                                  #use_DEG = config$analysis$use_DEG
                                 #)
#message(sprintf("%d common genes detected!", length(common_genes)))


#message("Starting meta-analysis...")

#combined_pData <- c(dna_pData, rna_pData)

#pData <- names(c(dna_matrices, norm_rna_mtxs))
#study <- names(c(dna_matrices, norm_rna_mtxs))
                                     
#meta_res <- generate_list_for_meta_analysis(
  #DNA = config$analysis$modalities$DNA,
  #RNA = config$analysis$modalities$RNA,
  #list_of_dna_mtx = dna_matrices,
  #list_of_rna_mtx = norm_rna_mtxs,
  #list_of_pData = combined_pData[pData],
  #study = study,
  #common_genes = common_genes
#)

#message("Meta-analysis completed!")
