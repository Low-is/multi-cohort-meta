library(yaml)
library(jsonlite)
library(COCONUT)
source("run-meta/functions/meta_analysis_functions.R")

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
message("Matrices loaded!")

rna_mtxs <- lapply(rna_matrices, function(x) x$expr)

lapply(rna_mtxs, function(x) str(x))

# Need to add code that filters matrices to match dimmensions of pData
all_matrices <- c(dna_matrices, rna_matrices)

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
