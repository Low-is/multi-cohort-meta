run_multicohort_meta <- function(
    DNA = FALSE,
    RNA = FALSE,
    list_of_dna_mtx = NULL,
    list_of_rna_mtx = NULL,
    list_of_dna_pheno = NULL,
    list_of_rna_pheno = NULL,
    gene_set = NULL,
    abs_ES_threshold = 0.02,
    fdr_threshold = 0.05,
    reference_level = "Control",
    het_p_threshold = NULL # <-- user can suppy NULL or numeric threshold
) {
  
  # --------------------------------------------------
  # 1. Input validation
  # --------------------------------------------------
  if (!DNA && !RNA) {
    stop("At least one of DNA or RNA must be TRUE.")
  }
  
  study_expr_list_dna  <- list()
  study_pheno_list_dna <- list()
  
  study_expr_list_rna  <- list()
  study_pheno_list_rna <- list()
  
  study_expr_list <- list()
  study_pheno_list <- list()
  
  # matrices_dna[c("GSE25504", "GSE69686", "GSE8586", "GSE32472")]
  # list(GSE138712 = sepsis_norm_gse138712_rna[["norm_counts"]],
  #      GSE125873 = bpd_norm_gse125873_rna[["norm_counts"]],
  #      GSE220135 = bpd_norm_gse220135_rna[["norm_counts"]])
  # 
  # sepsis_pData_dna[c("GSE25504", "GSE69686", "GSE8586", "GSE32472")]
  # sepsis_pData_rna[c("GSE138712", "GSE125873", "GSE220135")]
  
  if (DNA) {
    if (is.null(list_of_dna_mtx) || is.null(list_of_dna_pheno)) {
      stop("DNA selected but matrices or phenotype list missing.")
    }
    study_expr_list_dna  <- c(study_expr_list_dna, list_of_dna_mtx)
    study_pheno_list_dna <- c(study_pheno_list_dna, list_of_dna_pheno)
    # study_expr_list_dna  <- c(study_expr_list_dna, matrices_dna[c("GSE25504", "GSE69686", "GSE8586", "GSE32472")])
    # study_pheno_list_dna <- c(study_pheno_list_dna, sepsis_pData_dna[c("GSE25504", "GSE69686", "GSE8586", "GSE32472")])
  }
  
  if (RNA) {
    if (is.null(list_of_rna_mtx) || is.null(list_of_rna_pheno)) {
      stop("RNA selected but matrices or phenotype list missing.")
    }
    study_expr_list_rna  <- c(study_expr_list_rna, list_of_rna_mtx)
    study_pheno_list_rna <- c(study_pheno_list_rna, list_of_rna_pheno)
    # study_expr_list_rna  <- c(study_expr_list_rna, list(GSE138712 = sepsis_norm_gse138712_rna[["norm_counts"]],
    #                                                     GSE125873 = bpd_norm_gse125873_rna[["norm_counts"]],
    #                                                     GSE220135 = bpd_norm_gse220135_rna[["norm_counts"]])
    # )
    # study_pheno_list_rna <- c(study_pheno_list_rna, sepsis_pData_rna[c("GSE138712", "GSE125873", "GSE220135")])
  }
  
  # --------------------------------------------------
  # 2. Determine common genes
  # --------------------------------------------------
  if (!is.null(gene_set)) {
    common_genes <- gene_set
  } else {
    common_genes <- Reduce(intersect, lapply(c(study_expr_list_dna, study_expr_list_rna), rownames))
  }
  if (length(common_genes) == 0)
    stop("No shared genes found across selected studies.")
  
  study_expr_list <- lapply(c(study_expr_list_dna, study_expr_list_rna), function(m) m[common_genes, , drop = FALSE])
  study_pheno_list <- c(study_pheno_list_dna, study_pheno_list_rna)
  
  # --------------------------------------------------
  # 3. Standardize conditions
  # --------------------------------------------------
  standardize_condition <- function(cond) {
    message("Ensure 'condition' is a factor and control is first level.")
    
    if (!is.factor(cond)) {
      cond <- factor(cond)
    }
    
    if (length(levels(cond)) < 2) {
      stop("cond needs 2 levels.")
    }
    
    original_levels <- levels(cond)
    
    # Rename levels
    new_levels <- c(
      "Control",
      paste0("Case_", original_levels[2])
    )
    
    levels(cond) <- new_levels
    
    return(cond)
  }
  
  standardized_pheno <- lapply(study_pheno_list, function(p) {
    standardize_condition(p$condition)
  })
  
  combined_cases <- unlist(standardized_pheno)
  
  all_levels <- unique(as.character(combined_cases))
  
  all_pheno <- factor(combined_cases, levels = all_levels)
  
  # --------------------------------------------------
  # 4. Combine expression matrices
  # --------------------------------------------------
  batch <- factor(rep(
    names(study_expr_list),
    sapply(study_expr_list, ncol)
  ))
  
  all_expr <- do.call(cbind, study_expr_list)
  
  # --------------------------------------------------
  # 5. Batch correction
  # --------------------------------------------------
  design <- model.matrix(~ all_pheno)
  
  all_expr_bc <- limma::removeBatchEffect(
    all_expr,
    batch = batch,
    design = design
  )
  
  study <- factor(
    rep(
      names(study_expr_list),                  # study names
      sapply(study_expr_list, ncol)            # number of samples per study
    )
  )
  
  # --------------------------------------------------
  # 6. Effect size computation
  # --------------------------------------------------
  res <- effect.sizes.mg(expr = all_expr_bc, class = all_pheno, study = study)
  res$keys <- rownames(res)
  
  # Z-based p-value from eta²
  # z <- res$eta2 / res$se_eta2
  # res$p_value <- 2 * (1 - pnorm(abs(z)))
  # res$padj <- p.adjust(res$p_value, method = "BH")
  
  # Use ANOVA p-values directly
  res$padj <- p.adjust(res$p_value, method = "BH")
  
  
  rownames(res) <- res$keys
  # res <- res[complete.cases(res), ]
  
  # --------------------------------------------------
  # 7. Apply thresholds
  # --------------------------------------------------
  significant_genes <- res[
    abs(res$eta2) > abs_ES_threshold &
      res$padj < fdr_threshold,
  ]
  
  # --------------------------------------------------
  # 8. Return results
  # --------------------------------------------------
  return(list(
    all_results = res,
    significant_genes = significant_genes,
    batch_corrected_matrix = all_expr_bc,
    phenotype = all_pheno
  ))
}


# Usage
meta_res <- run_multicohort_meta(
  DNA = TRUE,
  RNA = TRUE,
  list_of_dna_mtx = matrices_dna[c("GSE25504", "GSE32472")],
  list_of_rna_mtx = list(GSE138712 = sepsis_norm_gse138712_rna[["norm_counts"]],
                         GSE220135 = bpd_norm_gse220135_rna[["norm_counts"]]),
  list_of_dna_pheno = sepsis_pData_dna[c("GSE25504", "GSE32472")],
  list_of_rna_pheno = sepsis_pData_rna[c("GSE138712", "GSE220135")],
  abs_ES_threshold = 0.02,
  fdr_threshold = 0.05
)

meta_res$significant_genes