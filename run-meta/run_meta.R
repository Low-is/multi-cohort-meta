library(yaml)
library(jsonlite)
library(COCONUT)
library(DESeq2)
source("run-meta/functions/meta_analysis_functions.R")
source("meta/functions/expr_mtx.R")

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
rna_matrices <- readRDS("meta/matrices/rna_matrices.rds") # the str is a list of 2: x$expr and x$pData
rna_matrices <- rna_matrices[!sapply(rna_matrices, is.null)] # Removing NULL single cell datasets
message("Matrices loaded!")

rna_pdata <- readRDS("meta/pdata/rna_pData_with_condition.rds")

lapply(rna_matrices, function(x) dim(x))
lapply(rna_pdata, function(x) dim(x))


#message("Getting norm RNA counts...")
#norm_rna_mtxs <- get_norm_RNA_counts(rna_matrices, pData = rna_pdata)
#message("Extracted norm RNA counts!")


# Need to add code that filters matrices to match dimmensions of pData

#all_matrices <- c(dna_matrices, rna_mtxs)
#lapply(all_matrices, function(x) head(rownames(x)))

# Find common genes across all studies being used for meta-analysis
#message("Searching for common genes...")
#common_genes <- Reduce(intersect, all_matrices)
#common_genes <- find_common_genes(DNA = config$analysis$modalities$DNA,
                                  #RNA = config$analysis$modalities$RNA,
                                  #list_of_dna_mtx = dna_mats,
                                  #list_of_rna_mtx = rna_mats,
                                  #use_DEG = config$analysis$use_DEG
                                 #)
#message(sprintf("%d common genes detected!", length(common_genes)))
