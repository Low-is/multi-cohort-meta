# Loading functions
source("meta/functions/get_pData.R")

# Load config
library(yaml)
config <- yaml::read_yaml("meta/config/config.yaml")

# Define search patterns
case_patterns <- config$case_patterns
control_patterns <- config$control_patterns

# Loading pData
dna_pData <- readRDS("meta/pdata/dna_pData.rds")
rna_pData <- readRDS("meta/pdata/rna_pData.rds")

# Apply condition detection + labeling
dna_pData_cond <- apply_condition_to_list(
  dna_pData,
  case_patterns,
  control_patterns
)

rna_pData_cond <- apply_condition_to_list(
  rna_pData,
  case_patterns,
  control_patterns
)

# Loading expression matrices
message("Loading expression matrices...")
dna_matrices <-  readRDS("meta/matrices/dna_matrices.rds")
dna_matrices <- dna_matrices[!sapply(dna_matrices, is.null)]

# Function to get rid of methylation experiements
is_meth <- function(genes) {
  mean(grepl("^cg", genes, ignore.case = TRUE)) > 0.05
}

dna_matrices <- dna_matrices[!sapply(dna_matrices, function(m) is_meth(rownames(m)))]

rna_matrices <- readRDS("meta/matrices/rna_matrices.rds") # the str is a list of 2: x$expr and x$pData
rna_matrices <- rna_matrices[!sapply(rna_matrices, is.null)] # Removing NULL single cell datasets
message("Matrices loaded!")

                                     
# Need to add code that filters matrices to match dimmensions of pData         
# Filtering list
dna_matrices <- dna_matrices[names(dna_pData_cond)] 
rna_matrices <- rna_matrices[names(rna_pData_cond)]

                                                                                           
message("Getting norm RNA counts...")
norm_rna_mtxs <- get_norm_RNA_counts(rna_matrices, pData = rna_pData_cond)
message("Extracted norm RNA counts!")

                                     
dna_matrices <- mapply(
  function(expr, pdata) {

    res <- resolve_matrix_names(
      expr_colnames = colnames(expr),
      pdata = pdata
    )

    expr <- expr[, res$pdata_rows, drop = FALSE]
    colnames(expr) <- res$matched_gsms

    expr
  },
  dna_matrices,
  dna_pData_cond,
  SIMPLIFY = FALSE
)

norm_rna_mtxs <- mapply(
  function(expr, pdata) {

    res <- resolve_matrix_names(
      expr_colnames = colnames(expr),
      pdata = pdata
    )

    expr <- expr[, res$pdata_rows, drop = FALSE]
    colnames(expr) <- res$matched_gsms

    expr
  },
  norm_rna_mtxs,
  rna_pData_cond,
  SIMPLIFY = FALSE
)


# Save results
saveRDS(dna_pData_cond, "meta/pdata/dna_pData_with_condition.rds")
saveRDS(rna_pData_cond, "meta/pdata/rna_pData_with_condition.rds")

message("✅ Condition column successfully added to all studies")

saveRDS(dna_matrices, "meta/matrices/dna_matrices.rds")
saveRDS(norm_rna_mtxs, "meta/matrices/norm_rna_mtxs.rds")

message("✅ Successfully saved updated matrices!")
