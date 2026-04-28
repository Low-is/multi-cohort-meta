############################################
# ---Meta-analysis function--- 
############################################

meta_results <- function(list_of_studies,
                         effect_thresh = 0.5,
                         fdr_thresh = 0.05,
                         het_thresh = 0.05) {
  
  # --------------------------
  # Compute effect sizes for each study
  # --------------------------
  list_of_es <- lapply(list_of_studies, effect.sizes)
  summary <- combine.effect.sizes(list_of_es)
  
  # --------------------------
  # Clean g and se.g
  # --------------------------
  bad_g    <- apply(is.na(summary$g)    | is.infinite(summary$g),    1, any)
  bad_se.g <- apply(is.na(summary$se.g) | is.infinite(summary$se.g), 1, any)
  keep <- !(bad_g | bad_se.g)
  
  summary$g    <- summary$g[keep, , drop = FALSE]
  summary$se.g <- summary$se.g[keep, , drop = FALSE]
  
  g_genes <- rownames(summary$g)
  
  # Clean pooled estimates to match filtered genes
  summary$pooled.estimates$genes <- rownames(summary$pooled.estimates)
  summary$pooled.estimates <- summary$pooled.estimates %>%
    dplyr::filter(genes %in% g_genes) %>%
    dplyr::filter(!apply(is.na(.) | is.infinite(as.matrix(.)), 1, any))
  
  # Avoid zero values in numeric columns
  num_cols <- sapply(summary$pooled.estimates, is.numeric)
  summary$pooled.estimates[num_cols] <- lapply(
    summary$pooled.estimates[num_cols],
    function(x) ifelse(x == 0, 1e-200, x)
  )
  
  rownames(summary$pooled.estimates) <- summary$pooled.estimates$genes
  
  # Extract main stats
  g     <- summary$g
  se.g  <- summary$se.g
  pool    <- summary$pooled.estimates$summary
  se.pool <- summary$pooled.estimates$se.summary
  p.het   <- summary$pooled.estimates$pval.het
  
  names(pool)    <- rownames(g)
  names(se.pool) <- rownames(g)
  names(p.het)   <- rownames(g)
  
  # Compute FDR
  common_genes <- Reduce(intersect, lapply(list_of_studies, function(s) rownames(s$expr)))
  study_pvals <- lapply(list_of_studies, function(study) {
    get.ttest.P(study$expr[common_genes, , drop=FALSE], study$class)[, "P.both"]
  })
  pval_matrix <- do.call(cbind, study_pvals)
  combined_pvals <- apply(pval_matrix, 1, function(pvec) metap::sumlog(pvec)$p)
  fdr <- p.adjust(combined_pvals, method = "fdr")
  
  summary$pooled.estimates$pooled_pval <- combined_pvals
  summary$pooled.estimates$FDR <- fdr
  
  # Confidence intervals
  lower_CI <- pool - 1.96 * se.pool
  upper_CI <- pool + 1.96 * se.pool
  
  # --------------------------
  # Robust gene selection using dynamic thresholds
  # --------------------------
  robust_genes <- character()
  for (gene in rownames(g)) {
    
    g_gene <- g[gene, ]
    
    strong_effect <- !is.na(pool[gene]) && abs(pool[gene]) > effect_thresh
    ci_consistent <- !is.na(lower_CI[gene]) && !is.na(upper_CI[gene]) &&
      ((lower_CI[gene] > 0) || (upper_CI[gene] < 0))
    hetero_ok     <- is.na(p.het[gene]) || p.het[gene] > het_thresh
    fdr_ok        <- !is.na(fdr[gene]) && fdr[gene] < fdr_thresh
    same_dir      <- all(!is.na(g_gene)) && (all(g_gene > 0) || all(g_gene < 0))
    
    if (strong_effect && ci_consistent && hetero_ok && fdr_ok && same_dir) {
      robust_genes <- c(robust_genes, gene)
    }
  }
  
  cat("Number of robust genes:", length(robust_genes), "\n")
  
  list(
    summary          = summary,
    pooled_estimates = summary$pooled.estimates,
    robust_genes     = robust_genes,
    g                = g,
    se.g             = se.g,
    pool             = pool,
    se.pool          = se.pool,
    p.het            = p.het,
    lower_CI         = lower_CI,
    upper_CI         = upper_CI,
    fdr              = fdr
  )
}

############################################
# ---Meta-analysis function--- 
############################################


############################################
# ---Meta-analysis function wrapper--- 
############################################

generate_list_for_meta_analysis <- function(
    DNA = FALSE, RNA = FALSE,
    list_of_dna_mtx = NULL,
    list_of_rna_mtx = NULL,
    list_of_pData = NULL,
    study = NULL,
    common_genes = NULL
) {
  
  list_of_studies <- list()
  
  # --------------------------
  # Detect common genes if not provided
  # --------------------------
  if (is.null(common_genes)) {
    message("➡️ No gene list provided. Automatically detecting common genes...")
    
    if (DNA && !RNA) {
      common_genes <- find_common_genes(DNA = TRUE, list_of_dna_mtx = list_of_dna_mtx)
    } else if (RNA && !DNA) {
      common_genes <- find_common_genes(RNA = TRUE, list_of_rna_mtx = list_of_rna_mtx)
    } else {
      common_genes <- find_common_genes(
        DNA = TRUE, list_of_dna_mtx = list_of_dna_mtx,
        RNA = TRUE, list_of_rna_mtx = list_of_rna_mtx
      )
    }
  }
  
  # --------------------------
  # Expression sanitization function
  # --------------------------
  sanitize_expr <- function(mat) {
    mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
    mat <- mat[apply(mat, 1, var) > 0, , drop = FALSE]
    mat
  }
  
  # --------------------------
  # Select the appropriate matrix list
  # --------------------------
  if (DNA && !RNA) {
    mtx_list <- list_of_dna_mtx
  } else if (RNA && !DNA) {
    mtx_list <- list_of_rna_mtx
  } else {
    # Combine DNA + RNA safely as list of lists
    mtx_list <- c(list_of_dna_mtx, list_of_rna_mtx)
  }
  
  # --------------------------
  # Loop over studies
  # --------------------------
  for (i in seq_along(study)) {
    mtx <- mtx_list[[i]]
    p   <- list_of_pData[[i]]
    
    genes_in_mtx <- intersect(common_genes, rownames(mtx))
    
    expr <- mtx[genes_in_mtx, , drop = FALSE]
    expr <- sanitize_expr(expr)
    
    if (nrow(expr) == 0) {
      warning(paste("⚠️ Study", study[i], "has no valid genes after filtering!"))
    }
    
    list_of_studies[[study[i]]] <- list(
      expr  = expr,
      pheno = p$condition,
      keys  = rownames(expr),
      class = as.numeric(ifelse(p$condition == levels(p$condition)[1], 0, 1))
    )
  }
  
  # --------------------------
  # Run meta-analysis
  # --------------------------
  meta <- meta_results(list_of_studies)
  
  # --------------------------
  # Return full structure
  # --------------------------
  structure(
    list(
      studies       = list_of_studies,
      robust_genes  = meta$robust_genes,
      pooled        = meta$pooled_estimates,
      meta          = meta
    ),
    class = "MetaAnalysis"
  )
}

############################################
# ---Meta-analysis function wrapper--- 
############################################



############################################
# ---Creating forest plots from meta-analysis--- 
############################################

plot_forest_interactive <- function(meta_obj, genes = NULL, xlim = c(-3, 3), colors = NULL) {
  if (is.null(colors)) colors <- list(
    box = "violetred",
    lines = "plum",
    summary = "mediumpurple",
    text = "black",
    axes = "black",
    zero = "black"
  )
  
  if (is.null(genes)) genes <- meta_obj$robust_genes
  if (length(genes) == 0) return(NULL)
  
  for (gene in genes) {
    g_gene <- meta_obj$g[gene, ]
    se.g_gene <- meta_obj$se.g[gene, ]
    pool <- meta_obj$pool[gene]
    se.pool <- meta_obj$se.pool[gene]
    study_names <- gsub("_g", "", names(g_gene))
    
    par(cex = 0.65)
    
    metaplot(
      g_gene, se.g_gene,
      labels = study_names,
      summn = pool,
      sumse = se.pool,
      sumnn = 1 / se.pool^2,
      summlabel = "Pooled Effect",
      xlab = "Standardized Mean Difference",
      ylab = "",
      xlim = xlim,
      main = bquote(italic(.(gene))),
      colors = meta.colors(
        box = colors$box,
        lines = colors$lines,
        summary = colors$summary,
        text = colors$text,
        axes = colors$axes,
        zero = colors$zero
      ),
      boxsize = 1,
      lty.random = 1,
      lwd.random = 2,
      zero = 0,
      col.zero = colors$zero,
      lty.zero = 3
    )
  }
  
  invisible(TRUE)
}


############################################
# ---Creating forest plots from meta-analysis--- 
############################################


############################################
# ---Export forest plots as PDF or PNG--- 
############################################

plot_forest_export <- function(meta_obj,
                               genes = NULL,
                               xlim = c(-3, 3),
                               colors = NULL,
                               out_dir = "forest_plots",
                               file_type = c("pdf", "png")) {
  
  file_type <- match.arg(file_type)
  if (is.null(colors)) colors <- list(
    box = "violetred",
    lines = "plum",
    summary = "mediumpurple",
    text = "black",
    axes = "black",
    zero = "black"
  )
  
  if (is.null(genes)) genes <- meta_obj$robust_genes
  if (length(genes) == 0) return(NULL)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  for (gene in genes) {
    file_path <- file.path(out_dir, paste0(gene, ".", file_type))
    
    # Open device
    if (file_type == "pdf") pdf(file_path, height = 3.5, width = 3.5)
    if (file_type == "png") png(file_path, height = 3.5, width = 3.5, units = "in", res = 300)
    
    # Plot
    g_gene <- meta_obj$g[gene, ]
    se.g_gene <- meta_obj$se.g[gene, ]
    pool <- meta_obj$pool[gene]
    se.pool <- meta_obj$se.pool[gene]
    study_names <- gsub("_g", "", names(g_gene))
    
    par(cex = 0.65)
    
    metaplot(
      g_gene, se.g_gene,
      labels = study_names,
      summn = pool,
      sumse = se.pool,
      sumnn = 1 / se.pool^2,
      summlabel = "Pooled Effect",
      xlab = "Standardized Mean Difference",
      ylab = "",
      xlim = xlim,
      main = bquote(italic(.(gene))),
      colors = meta.colors(
        box = colors$box,
        lines = colors$lines,
        summary = colors$summary,
        text = colors$text,
        axes = colors$axes,
        zero = colors$zero
      ),
      boxsize = 1,
      lty.random = 1,
      lwd.random = 2,
      zero = 0,
      col.zero = colors$zero,
      lty.zero = 3
    )
    
    dev.off()
  }
  
  invisible(TRUE)
}


############################################
# ---Export forest plots as PDF or PNG--- 
############################################