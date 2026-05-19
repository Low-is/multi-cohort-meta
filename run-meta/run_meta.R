library(yaml)
library(jsonlite)
library(COCONUT)
library(DESeq2)
library(data.table)
library(dplyr)
library(rmeta)
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


# ----------------------------
# NORMALIZATION (GLOBAL SAFE)
# ----------------------------
normalize <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9 ]", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}


# Need to add code that filters matrices to match dimmensions of pData         
# Filtering list
dna_matrices <- dna_matrices[!is.na(names(dna_matrices))]
dna_pData <- dna_pData[names(dna_matrices)]

dna_matrices <- mapply(function(mat, pd) {
  if (is.null(mat) || is.null(pd)) return(NULL)

  common_gsm <- intersect(colnames(mat), pd$gsm)
  mat <- mat[, common_gsm, drop = FALSE]

  return(mat)
},
                       dna_matrices, dna_pData, SIMPLIFY = FALSE)

                                     

dna_matrices <- lapply(dna_matrices, function(x) {
  if (is.null(x)) return(NULL)
  x[complete.cases(x), , drop = FALSE]
})

valid_studies <- sapply(dna_pData, function(pd) {
  tab <- table(pd$condition)

  length(tab) >= 2 %% all(tab >= 2)
})

valid_studies[is.na(valid_studies)] <- FALSE


dna_matrices <- dna_matrices[valid_studies]
dna_pData <- dna_pData[valid_studies]

names(c(dna_matrices, rna_matrices))
                                                                                                 
# Find common genes across all studies being used for meta-analysis
message("Searching for common genes...")
common_genes <- find_common_genes(DNA = config$analysis$modalities$DNA,
                                  RNA = config$analysis$modalities$RNA,
                                  list_of_dna_mtx = dna_matrices,
                                  list_of_rna_mtx = rna_matrices,
                                  use_DEG = config$analysis$use_DEG
                                 )
message(sprintf("%d common genes detected!", length(common_genes)))


message("Starting meta-analysis...")

meta_list <- list(
  GSE8586 = list(
    expr = dna_matrices[["GSE8586"]][common_genes, ],  
    pheno = dna_pData[["GSE8586"]]$condition, 
    keys = rownames(dna_matrices[["GSE8586"]][common_genes, ]),  
    class = factor(as.numeric(ifelse(dna_pData[["GSE8586"]]$condition == "Control", 0, 1), levels = c(0, 1)))
  ),
  GSE32472 = list(
    expr = dna_matrices[["GSE32472"]][common_genes, ],  
    pheno = dna_pData[["GSE32472"]]$condition, 
    keys = rownames(dna_matrices[["GSE32472"]][common_genes, ]),  
    class = factor(as.numeric(ifelse(dna_pData[["GSE32472"]]$condition == "Control", 0, 1), levels = c(0, 1)))
  ), 
  GSE188944 = list(
    expr = dna_matrices[["GSE188944"]][common_genes, ],  
    pheno = dna_pData[["GSE188944"]]$condition, 
    keys = rownames(dna_matrices[["GSE188944"]][common_genes, ]),  
    class = factor(as.numeric(ifelse(dna_pData[["GSE188944"]]$condition == "Control", 0, 1), levels = c(0, 1)))
  ),
  GSE225881 = list(
    expr = dna_matrices[["GSE225881"]][common_genes, ],  
    pheno = dna_pData[["GSE225881"]]$condition, 
    keys = rownames(dna_matrices[["GSE225881"]][common_genes, ]),  
    class = factor(as.numeric(ifelse(dna_pData[["GSE225881"]]$condition == "Control", 0, 1), levels = c(0, 1)))
  ),
  GSE108754 = list(
    expr = dna_matrices[["GSE108754"]][common_genes, ],  
    pheno = dna_pData[["GSE108754"]]$condition, 
    keys = rownames(dna_matrices[["GSE108754"]][common_genes, ]),  
    class = factor(as.numeric(ifelse(dna_pData[["GSE108754"]]$condition == "Control", 0, 1), levels = c(0, 1)))
  ),
  GSE108756 = list(
    expr = dna_matrices[["GSE108756"]][common_genes, ],  
    pheno = dna_pData[["GSE108756"]]$condition, 
    keys = rownames(dna_matrices[["GSE108756"]][common_genes, ]),  
    class = factor(as.numeric(ifelse(dna_pData[["GSE108756"]]$condition == "Control", 0, 1), levels = c(0, 1)))
  ),
  GSE219156 = list(
    expr = rna_matrices[["GSE219156"]][common_genes, ],  
    pheno = rna_pData[["GSE219156"]]$condition, 
    keys = rownames(rna_matrices[["GSE219156"]][common_genes, ]),  
    class = factor(as.numeric(ifelse(rna_pData[["GSE219156"]]$condition == "Control", 0, 1), levels = c(0, 1)))
  ),
  GSE106910 = list(
    expr = rna_matrices[["GSE106910"]][common_genes, ],  
    pheno = rna_pData[["GSE106910"]]$condition, 
    keys = rownames(rna_matrices[["GSE106910"]][common_genes, ]),  
    class = factor(as.numeric(ifelse(rna_pData[["GSE106910"]]$condition == "Control", 0, 1), levels = c(0, 1)))
  ),
  GSE125873 = list(
    expr = rna_matrices[["GSE125873"]][common_genes, ],  
    pheno = rna_pData[["GSE125873"]]$condition, 
    keys = rownames(rna_matrices[["GSE125873"]][common_genes, ]),  
    class = factor(as.numeric(ifelse(rna_pData[["GSE125873"]]$condition == "Control", 0, 1), levels = c(0, 1)))
  ),
  GSE220135 = list(
    expr = rna_matrices[["GSE220135"]][common_genes, ],  
    pheno = rna_pData[["GSE220135"]]$condition, 
    keys = rownames(rna_matrices[["GSE220135"]][common_genes, ]),  
    class = factor(as.numeric(ifelse(rna_pData[["GSE220135"]]$condition == "Control", 0, 1), levels = c(0, 1))))
)

gse8586_es <- effect.sizes(meta_list$GSE8586)
gse32472_es <- effect.sizes(meta_list$GSE32472)
gse188944_es <- effect.sizes(meta_list$GSE188944)
gse225881_es <- effect.sizes(meta_list$GSE225881)
gse108754_es <- effect.sizes(meta_list$GSE108754)
gse108756_es <- effect.sizes(meta_list$GSE108756)
gse219156_es <- effect.sizes(meta_list$GSE219156)                                    
gse106910_es <- effect.sizes(meta_list$GSE106910)
gse125873_es <- effect.sizes(meta_list$GSE125873)
gse220135_es <- effect.sizes(meta_list$GSE220135)

                                     
list.of.effects <- list(GSE8586 = gse8586_es,
                        GSE32472 = gse32472_es,
                        GSE188944 = gse188944_es,
                        GSE225881 = gse225881_es,
                        GSE108754 = gse108754_es,
                        GSE108756 = gse108756_es,
                        GSE219156 = gse219156_es,
                        GSE106910 = gse106910_es,
                        GSE125873 = gse125873_es,
                        GSE220135 = gse220135_es)

summary <- combine.effect.sizes(list.of.effects)  

summary$g <- summary$g[!apply(is.na(summary$g) | is.infinite(summary$g), 1, any), ]
summary$se.g <- summary$se.g[!apply(is.na(summary$se.g) | is.infinite(summary$se.g), 1, any), ]
g_genes <- rownames(summary$g)
    
summary$pooled.estimates$genes <- rownames(summary$pooled.estimates)
summary$pooled.estimates <- summary$pooled.estimates %>%
                                     dplyr::filter(genes %in% g_genes)
    
summary$pooled.estimates <- summary$pooled.estimates %>%
                                     dplyr::filter(!apply(is.na(.) | is.infinite(as.matrix(.)), 1, any))
    
summary$pooled.estimates <- as.data.frame(lapply(summary$pooled.estimates, function(x) ifelse(x == 0, 1e-200, x)))
rownames(summary$pooled.estimates) <- g_genes

# Extracting individual estimates                                                 
g <- summary$g
se.g <- summary$se.g

# Extracting pooled estimates
pool    <- summary$pooled.estimates[, "summary"]
names(pool) <- rownames(g)
    
se.pool <- summary$pooled.estimates[, "se.summary"]
names(se.pool) <- rownames(g)

 # Label for forestplots                                               
x.label <- "Standardized Mean Difference (log2 scale)"
    
# Adding FDR corrected p-values
summary$pooled.estimates <- summary$pooled.estimates %>%
                                                 dplyr::mutate(
                                                   FDR = p.adjust(p.value, method = "BH")
                                                 )
                                                     
summary$pooled.estimates
#combined_pData <- c(dna_pData, rna_pData)

#pData <- names(c(dna_matrices, rna_matrices))
#study <- names(c(dna_matrices, rna_matrices))
                                     
#meta_res <- generate_list_for_meta_analysis(
  #DNA = config$analysis$modalities$DNA,
  #RNA = config$analysis$modalities$RNA,
  #list_of_dna_mtx = dna_matrices,
  #list_of_rna_mtx = rna_matrices,
  #list_of_pData = combined_pData,
  #study = study,
  #common_genes = common_genes
#)

#message("Meta-analysis completed!")
