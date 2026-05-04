############################################
# ---Finding common genes across studies---
############################################

find_common_genes <- function(
    DNA = FALSE, RNA = FALSE,
    list_of_dna_mtx = NULL,
    list_of_rna_mtx = NULL,
    use_DEG = FALSE,
    dna_deg_list = NULL,
    rna_deg_list = NULL,
    sig_label = "Significant"
) {
  
  # --- Helper: intersect genes across a list of objects ---
  intersection_from_lists <- function(lst) {
    Reduce(intersect, lst)
  }
  
  # ============================================================
  #       MODE 1: Use ONLY rownames (original behavior)
  # ============================================================
  if (!use_DEG) {
    
    if (DNA && !RNA) {
      if (!is.list(list_of_dna_mtx))
        stop("DNA matrices need to be in a list")
      
      return(intersection_from_lists(lapply(list_of_dna_mtx, rownames)))
    }
    
    if (RNA && !DNA) {
      if (!is.list(list_of_rna_mtx))
        stop("RNA matrices need to be in a list")
      
      return(intersection_from_lists(lapply(list_of_rna_mtx, rownames)))
    }
    
    if (DNA && RNA) {
      if (!is.list(list_of_dna_mtx) || !is.list(list_of_rna_mtx))
        stop("Both list_of_dna_mtx and list_of_rna_mtx must be lists of matrices.")
      
      list_of_mtx <- c(list_of_dna_mtx, list_of_rna_mtx)
      return(intersection_from_lists(lapply(list_of_mtx, rownames)))
    }
    
    return(NULL)
  }
  
  # ============================================================
  #       MODE 2: Use DE results ("Significance" == sig_label)
  # ============================================================
  if (use_DEG) {
    
    # --- DNA DEGs ---
    dna_sig_genes <- NULL
    if (DNA) {
      if (!is.list(dna_deg_list))
        stop("dna_deg_list must be a list of DEG data frames.")
      
      dna_sig_genes <- lapply(dna_deg_list, function(df) {
        if (!"Significance" %in% colnames(df))
          stop("DNA DEG result missing a 'Significance' column.")
        
        rownames(df)[df$Significance == sig_label]
      })
      dna_sig_genes <- intersection_from_lists(dna_sig_genes)
    }
    
    # --- RNA DEGs ---
    rna_sig_genes <- NULL
    if (RNA) {
      if (!is.list(rna_deg_list))
        stop("rna_deg_list must be a list of DEG data frames.")
      
      rna_sig_genes <- lapply(rna_deg_list, function(df) {
        if (!"Significance" %in% colnames(df))
          stop("RNA DEG result missing a 'Significance' column.")
        
        rownames(df)[df$Significance == sig_label]
      })
      rna_sig_genes <- intersection_from_lists(rna_sig_genes)
    }
    
    # Combine DNA + RNA intersections if both requested
    if (DNA && RNA) {
      return(intersect(dna_sig_genes, rna_sig_genes))
    }
    
    if (DNA) return(dna_sig_genes)
    if (RNA) return(rna_sig_genes)
  }
  
  return(NULL)
}

############################################
# ---Finding common genes across studies---
############################################


############################################
# ---Meta-analysis function--- 
############################################

meta_results <- function(list_of_studies) {
    
    list_of_es <- lapply(list_of_studies, function(study) {
      effect.sizes(study)
    })
    
    summary <- combine.effect.sizes(list_of_es)
    
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
    
    g <- summary$g
    se.g <- summary$se.g
    
    pool    <- summary$pooled.estimates[, "summary"]
    names(pool) <- rownames(g)
    
    se.pool <- summary$pooled.estimates[, "se.summary"]
    names(se.pool) <- rownames(g)
    
    x.label <- "Standardized Mean Difference (log2 scale)"
    
    # Adding FDR corrected p-values
    summary$pooled.estimates <- summary$pooled.estimates %>%
      dplyr::mutate(
        FDR = p.adjust(p.value, method = "BH")
      )
    

    
    consistent_genes <- c()
    
    for (gene in rownames(g)) {
      
      g_gene  <- g[gene, ]
      se_gene <- se.g[gene, ]
      
      pooled   <- pool[gene]
      se_pooled <- se.pool[gene]
      
      # 95% CI
      lower_CI <- pooled - 1.96 * se_pooled
      upper_CI <- pooled + 1.96 * se_pooled
      
      # Criterion 1: pooled effect must be ≥ |0.5|
      strong_effect <- abs(pooled) >= 0.5
      
      # Criterion 2: CI does NOT cross zero
      ci_significant <- (lower_CI > 0) || (upper_CI < 0)
      
      # Criterion 3 (optional): All study effects trend same direction
      # (keep this if you want it, or remove if not needed)
      consistent_direction <- (all(g_gene > 0)) || (all(g_gene < 0))
      
      # Final rule
      if (strong_effect && ci_significant && consistent_direction) {
        consistent_genes <- c(consistent_genes, gene)
      }
    }

    fp_data <- summary$pooled_estimates[consistent_genes, ] %>%
      dplyr::mutate(
          logOR = summary,
          SE = se.summary,
          OR = exp(logOR),
          lower = exp(logOR - 1.96 * SE),
          upper = exp(logOR + 1.96 * SE),
          FDR = round(FDR, 3)
      )

    fp_data <- fp_data %>%
      dplyr::mutate(
          Gene = consistent_genes,
          CI = sprintf("%.2f - %.2f", lower, upper),
          OR_txt = sprintf("%.2f", OR)
      ) %>%
      dplyr::arrange(desc(OR), Gene)

    label_mat <- cbind(
        as.character(fp_data$Gene),
        fp_data$OR_txt,
        fp_data$CI,
        fp_data$FDR
    )
    colnames(label_mat) <- c("Gene", "OR", "CI", "FDR")
    png("run-meta/output/fp.png",
       width = 10,
       height = 14,
       units = "in",
       res = 800)
    forestplot(
        labeltext = label_mat,
        mean = fp_data$OR,
        lower = fp_data$lower,
        upper = fp_data$upper,

        boxsize = 0.4,
        clip = c(0.1, 10),
        xticks = c(0.1, 0.5, 1, 2, 5, 10),

        xlog = TRUE,
        zero = 1,

        graphwidth = unit(0.5, "npc"),

        colgap = unit(6, "mm"),
        lineheight = unit(6, "mm"),

        col = fpColors(
            box = "royalblue",
            line = "darkblue",
            summary = "royalblue"
        ),
        txt_gp = fpTxtGp(
            ticks = gpar(cex = 1.3),
            xlab = gpar(cex = 1.5)
        )
    ) |>
        fp_add_header(
            Gene = c("", "Gene"),
            OR = c("", "OR"),
            CI = c("", "95% CI"),
            FDR = c("", "FDR")
        ) |>
        fp_set_zebra_style("#EFEFEF") |>
        fp_set_favors(low = "Likely No BPD",
                     high = "Likely BPD",
                     txt_gp = gapr(cex = 1.5),
                     arrows = FALSE)
    dev.off()
    
    # Flatten results for easier access
    meta_flat <- list(
      summary          = summary,
      pooled_estimates = summary$pooled.estimates,
      consistent_genes = consistent_genes
    )
    
    return(meta_flat)
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


  sanitize_expr <- function(mat) {
    mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
    mat <- mat[apply(mat, 1, var) > 0, , drop = FALSE]
    mat
  }


  for (i in seq_along(study)) {
    mtx <- c(list_of_dna_mtx, list_of_rna_mtx)[[i]]
    p   <- list_of_pData[[i]]

    genes_in_mtx <- intersect(common_genes, rownames(mtx))

    expr <- mtx[genes_in_mtx, , drop = FALSE]
    expr <- sanitize_expr(expr)

    list_of_studies[[study[i]]] <- list(
      expr  = expr,
      pheno = p$condition,
      keys  = rownames(expr),
      class = as.numeric(ifelse(p$condition == levels(p$condition)[1], 0, 1))
    )

  }

  meta <- meta_results(list_of_studies)

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
# ---Meta-analysis internal functions--- 
############################################

# -------------------------------------------
# Performing t-test 
# ------------------------------------------- 
get.ttest.P <- function(mat, g) {
  # mat: genes x samples
  # g: numeric vector of class labels (0/1 or 1/2)
  
  # Compute t-statistics with mt.teststat
  tstat <- mt.teststat(mat, g, test = "t") 
  
  # degrees of freedom
  df <- length(g) - 2
  
  # Compute p-values from t-statistics
  P.both <- 2 * pt(abs(tstat), df = df, lower.tail = FALSE)
  P.down <- pt(tstat, df = df, lower.tail = TRUE)
  P.up   <- pt(tstat, df = df, lower.tail = FALSE)
  
  out <- cbind(P.both, P.down, P.up)
  rownames(out) <- rownames(mat)
  return(out)
}
# -------------------------------------------
# Performing t-test 
# ------------------------------------------- 


# -------------------------------------------
# Pooling estimates - inverse variance weighting
# -------------------------------------------
pool.inverseVar <- function( g, se.g, method ){
  stopifnot( identical( rownames(g), rownames(se.g) ) )
  out <- matrix( nr=nrow(g), nc=8,
                 dimnames=list( rownames(g), c("n.studies", "summary", "se.summary", "tau2", "p.value", "Q", "df", "pval.het") ) )
  
  for(j in 1:nrow(g)){
    
    e  <- cleanNA(    g[j, ] )
    se <- cleanNA( se.g[j, ] )
    n  <- length(e)
    
    if(n==1){
      summ <- e;   se.summ <- se;   tau2 <- NA
      Q.het = NA
      df.het = NA
      pval.het = NA
    } else {
      fit <- meta.summaries(e, se, method = method)
      summ <- fit$summary
      se.summ <- fit$se.summary
      tau2 <- ifelse( method=="fixed", NA, fit$tau2 )
      Q.het = fit$het[1]
      df.het = fit$het[2]
      pval.het = fit$het[3]
      rm(fit)
    }
    
    pval     <- 2*pnorm( abs(summ/se.summ), lower.tail=FALSE )
    
    out[j, ] <- c(n, summ, se.summ, tau2, pval, Q.het, df.het, pval.het)
    rm(e, se, n, summ, se.summ, tau2, pval, Q.het, df.het, pval.het)
  }
  return(out)
}

# -------------------------------------------
# Pooling estimates - inverse variance weighting
# -------------------------------------------

############################################
# ---Meta-analysis internal functions--- 
############################################



############################################
# ---Leave One Out Cross-validation--- 
############################################

leaveOneOutMetaAnalysis <- function(studies, fdr_thresh, min_summary, hetero_pval_thresh) {
  loo_results <- vector("list", length(studies))
  names(loo_results) <- names(studies)
  
  for (i in seq_along(studies)) {
    cat("➡️ Leave-one-out iteration: leaving out ", names(studies)[i], "...\n", sep="")
    subset_studies <- studies[-i]
    
    # Compute study-level effect sizes
    all_ES <- lapply(subset_studies, effect.sizes)
    
    # Combine into pooled estimates
    output_RE <- combine.effect.sizes(all_ES)
    pe <- output_RE$pooled.estimates
    
    # Adjust FDRs
    pe$FDR_effectsize <- p.adjust(pe$p.value, method = "fdr")
    
    # Fisher combination of p-values
    all_Pvals <- lapply(subset_studies, ttest.Pvalues)
    fisher_out <- adjust.fisher(sum.of.logs(all_Pvals))
    
    # Ensure it is a data.frame
    if (is.vector(fisher_out) || is.atomic(fisher_out)) {
      # turn named vector into 1‑column data.frame
      fisher_out <- data.frame(F.Qval.up = as.numeric(fisher_out),
                               stringsAsFactors = FALSE)
      if (!is.null(names(fisher_out$F.Qval.up))) {
        rownames(fisher_out) <- names(fisher_out$F.Qval.up)
      }
    }
    
    # Now it is safe to subset by rownames
    fisher_out <- fisher_out[rownames(pe), , drop = FALSE]
    
    pe$FDR_Fisher <- fisher_out$F.Qval.up
    
    # Select genes for this LOO iteration
    selected_genes <- rownames(pe)[
      (pe$FDR_effectsize < fdr_thresh | pe$FDR_Fisher < fdr_thresh) &
        (pe$pval.het > hetero_pval_thresh) &
        (abs(pe$summary) > min_summary)
    ]
    
    loo_results[[i]] <- list(
      left_out = names(studies)[i],
      pooled_estimates = pe,
      selected_genes = selected_genes
    )
    
    cat("✅ Selected ", length(selected_genes), " genes\n", sep="")
  }
  
  # Robust genes = intersection across all LOO iterations
  robust_genes <- Reduce(intersect, lapply(loo_results, `[[`, "selected_genes"))
  cat("🎯 Total robust genes: ", length(robust_genes), "\n", sep="")
  
  list(
    loo_results = loo_results,
    robust_genes = robust_genes
  )
}

############################################
# ---Leave One Out Cross-validation--- 
############################################


############################################
# ---Greedy Forward Search--- 
############################################

# 1️⃣ Convert meta-analysis results into discovery.genes list for forward search
convertDiscoveryListToGenes_noDirection <- function(studies, genes) {
  discovery.genes <- lapply(studies, function(study) {
    expr <- study$expr[genes, , drop = FALSE]
    
    # Shift to non-negative
    xmin <- min(expr, na.rm = TRUE)
    if (xmin < 0) expr <- expr + abs(xmin)
    
    list(
      genes = expr,
      class = study$class
    )
  })
  return(discovery.genes)
}


# 2️⃣ Forward search function
forwardSearchGreedy <- function(
    discovery.genes,
    candidate_genes,
    forwardThresh = 0.01
) {
  selected_genes <- character()
  auc.current <- 0
  remaining <- candidate_genes
  
  repeat {
    res <- lapply(remaining, function(g) {
      genes_try <- c(selected_genes, g)
      auc_new <- getWeightedAUCsGenesList(discovery.genes, genes_try)
      data.frame(
        gene = g,
        auc = auc_new,
        diff = auc_new - auc.current
      )
    })
    res <- do.call(rbind, res)
    best <- max(res$diff, na.rm = TRUE)
    cat("Next best ΔAUC:", best, "\n")
    if (best < forwardThresh) break
    best_row <- res[res$diff == best, ][1, ]
    selected_genes <- c(selected_genes, best_row$gene)
    remaining <- setdiff(remaining, best_row$gene)
    auc.current <- best_row$auc
    cat("Adding", best_row$gene, "→ AUC =", auc.current, "\n")
  }
  
  list(
    selected_genes = selected_genes,
    final_auc = auc.current
  )
}

# 3️⃣ Weighted AUC helper
getWeightedAUCsGenesList <- function(discovery.genes, genes) {
  aucs <- sapply(discovery.genes, function(GEM) {
    X <- GEM$genes[genes, , drop = FALSE]
    xmin <- min(X, na.rm = TRUE)
    if (xmin < 0) X <- X + abs(xmin)
    score <- colMeans(X)
    efficientAUC(GEM$class, score)
  })
  weights <- sapply(discovery.genes, function(GEM) length(GEM$class))
  sum(aucs * weights) / sum(weights)
}

# 4️⃣ Efficient AUC function
efficientAUC <- function(labels, predictions) {
  levels <- sort(unique(labels))
  labels <- ordered(labels, levels = levels)
  n.pos <- sum(labels == levels[2])
  n.neg <- sum(labels == levels[1])
  pred.order <- order(predictions, decreasing = TRUE)
  predictions.sorted <- predictions[pred.order]
  tp <- cumsum(labels[pred.order] == levels[2])
  fp <- cumsum(labels[pred.order] == levels[1])
  dups <- rev(duplicated(rev(predictions.sorted)))
  tp <- c(0, tp[!dups])
  fp <- c(0, fp[!dups])
  fn <- n.pos - tp
  tn <- n.neg - fp
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  roc <- data.frame(x = fpr, y = tpr)
  roc <- roc[order(roc$x, roc$y),]
  i <- 2:nrow(roc)
  auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
  return(auc)
}

############################################
# ---Greedy Forward Search--- 
############################################


############################################
# ---Merging Expression Matrices--- 
############################################

merge_and_batch_correct_expression <-  function(dna_list,
                                                rna_list,
                                                dna_pData,
                                                rna_pData) {
  
  # Safety checks
  if (!is.list(dna_list) || !is.list(rna_list)) {
    stop("dna_list and rna_list must be named lists of matrices")
  }
  if (!is.list(dna_pData) || !is.list(rna_pData)) {
    stop("dna_pData and rna_pData must be named lists of data frames")
  }
  
  # ---- Find common genes ----
  common_genes <- find_common_genes(
    DNA = TRUE, RNA = TRUE,
    list_of_dna_mtx = dna_list,
    list_of_rna_mtx = rna_list
  )
  
  # ---- Subset matrices ----
  dna_list_sub <- lapply(dna_list, function(m) m[common_genes, , drop = FALSE])
  rna_list_sub <- lapply(rna_list, function(m) m[common_genes, , drop = FALSE])
  
  # ---- Merge matrices ----
  all_mtx <- c(dna_list_sub, rna_list_sub)
  merged_expr <- do.call(cbind, all_mtx)
  
  # ---- Build batch vector ----
  dna_batches <- rep(seq_along(dna_list_sub),
                     times = sapply(dna_list_sub, ncol))
  
  rna_offset <- length(dna_list_sub)
  rna_batches <- rep(seq_len(length(rna_list_sub)) + rna_offset,
                     times = sapply(rna_list_sub, ncol))
  
  batch <- c(dna_batches, rna_batches)
  
  # ---- Extract condition vector from pData ----
  dna_conditions <- unlist(lapply(dna_pData, function(p) {
    if (!"condition" %in% colnames(p)) {
      stop("Missing 'condition' column in DNA pData")
    }
    as.character(p$condition)
  }))
  
  rna_conditions <- unlist(lapply(rna_pData, function(p) {
    if (!"condition" %in% colnames(p)) {
      stop("Missing 'condition' column in RNA pData")
    }
    as.character(p$condition)
  }))
  
  condition <- factor(c(dna_conditions, rna_conditions))
  
  # ---- Design matrix to protect biology ----
  design <- model.matrix(~ condition)
  
  # ---- Batch correction ----
  corrected_expr <- limma::removeBatchEffect(
    merged_expr,
    batch = batch,
    design = design
  )
  
  # ---- Z-score scaling per gene ----
  scaled_expr <- t(scale(t(corrected_expr)))
  
  return(list(
    expr_matrix = scaled_expr,
    batch = batch,
    condition = condition,
    genes = common_genes
  ))
}

############################################
# ---Merging Expression Matrices--- 
############################################ 


############################################
# ---Wilcoxon--- 
############################################ 

# ### Plots wilcox of combined matrices ####
# plot_wilcox_genes <- function(combined_mtx_res,
#                               genes_of_interest,
#                               save_plot = FALSE,
#                               filename  = "wilcox_boxplots.png",
#                               width     = 8,
#                               height    = 6,
#                               dpi       = 600) {
# 
#     # Extract components
#     expr_mat  <- combined_mtx_res$expr_matrix
#     condition <- combined_mtx_res$condition
# 
#     # Keep only valid genes
#     genes_of_interest <- intersect(genes_of_interest, rownames(expr_mat))
#     if (length(genes_of_interest) == 0) {
#       stop("No valid genes found.")
#     }
# 
#     # Subset matrix
#     sub_expr <- expr_mat[genes_of_interest, , drop = FALSE]
# 
#     # Long-format dataframe
#     df <- data.frame(
#       gene       = rep(rownames(sub_expr), each = ncol(sub_expr)),
#       expression = as.vector(as.matrix(sub_expr)),
#       condition  = rep(condition, times = nrow(sub_expr))
#     )
# 
#     # Compute Wilcoxon p-values per gene
#     pvals <- df %>%
#       dplyr::group_by(gene) %>%
#       dplyr::summarise(
#         p_value = wilcox.test(expression ~ condition)$p.value,
#         .groups = "drop"
#       ) %>%
#       dplyr::mutate(
#         sig = dplyr::case_when(
#           p_value <= 0.001 ~ "***",
#           p_value <= 0.01  ~ "**",
#           p_value <= 0.05  ~ "*",
#           TRUE             ~ "ns"
#         ),
#         label = paste0("p = ", signif(p_value, 4), " ", sig)
#       )
# 
#     # Build annotation dataframe like your example
#     y_max <- df %>%
#       dplyr::group_by(gene) %>%
#       dplyr::summarise(y = max(expression, na.rm = TRUE) * 1.15)
# 
#     custom_p_L <- dplyr::left_join(pvals, y_max, by = "gene") %>%
#       dplyr::rename(p = label)
# 
#     # Create plot
#     p <- ggplot(df, aes(x = condition, y = expression, fill = condition)) +
#       geom_boxplot(outlier.shape = NA, alpha = 0.7) +
#       geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
#       facet_wrap(~ gene, scales = "free_y") +
#       geom_text(
#         data = custom_p_L,
#         aes(x = 1.5, y = y, label = p),
#         inherit.aes = FALSE,
#         size = 3.5
#       ) +
#       labs(x = NULL, y = "Scaled Expression") +
#       theme_bw(base_size = 12) +
#       theme(panel.grid = element_blank()) +
#       theme(
#         legend.position  = "none",
#         strip.background = element_blank(),
#         strip.text       = element_text(face = "bold")
#       )
# 
#     # Save if requested
#     if (save_plot) {
#       ggsave(
#         filename = filename,
#         plot     = p,
#         dpi      = dpi,
#         width    = width,
#         height   = height,
#         units    = "in"
#       )
#     }
# 
#     # Return both plot and p-values table
#     return(list(
#       plot = p,
#       p_values = custom_p_L
#     ))
#   }
# ### Plots wilcox of combined matrices ####

############################################
# ---Wilcoxon--- 
############################################ 


############################################
# ---T-test--- 
############################################ 

# plot_ttest_genes <- function(list_of_mtx,
#                                list_of_pData,
#                                genes_of_interest = NULL,
#                                save_plots = FALSE,
#                                prefix = "dataset_boxplot",
#                                width = 8,
#                                height = 6,
#                                dpi = 600) {
#     
#     results_list <- list()
#     
#     # Loop through each dataset separately
#     for (i in seq_along(list_of_mtx)) {
#       
#       mtx <- list_of_mtx[[i]]
#       pdata <- list_of_pData[[i]]
#       dataset_name <- names(list_of_mtx)[i]
#       if (is.null(dataset_name)) dataset_name <- paste0("Dataset_", i)
#       
#       # Ensure condition factor
#       condition <- factor(pdata$condition)
#       if (nlevels(condition) != 2)
#         stop(paste0("Dataset ", dataset_name, " does not have exactly 2 condition groups."))
#       
#       # ------------------------------
#       # 1. Z-score normalization
#       # ------------------------------
#       # mtx_z <- t(scale(t(mtx)))  # z-score per gene
#       
#       # ------------------------------
#       # 2. Select genes
#       # ------------------------------
#       if (!is.null(genes_of_interest)) {
#         genes_keep <- intersect(genes_of_interest, rownames(mtx))
#         if (length(genes_keep) == 0)
#           stop(paste0("No requested genes found in ", dataset_name))
#       } else {
#         genes_keep <- rownames(mtx)  # use all genes
#       }
#       
#       mtx_sub <- mtx[genes_keep, , drop = FALSE]
#       
#       # ------------------------------
#       # 3. Build long dataframe
#       # ------------------------------
#       df <- data.frame(
#         gene       = rep(rownames(mtx_sub), each = ncol(mtx_sub)),
#         expression = as.vector(as.matrix(mtx_sub)),
#         condition  = rep(condition, times = nrow(mtx_sub))
#       )
#       
#       # ------------------------------
#       # 4. Compute t-test p-values per gene
#       # ------------------------------
#       pvals <- df %>%
#         dplyr::group_by(gene) %>%
#         dplyr::summarise(
#           p_value = t.test(expression ~ condition, var.equal = FALSE)$p.value, # changed 'var.equal = TRUE' to FALSE
#           .groups = "drop"
#         ) %>%
#         dplyr::mutate(
#           sig = dplyr::case_when(
#             p_value <= 0.001 ~ "***",
#             p_value <= 0.01  ~ "**",
#             p_value <= 0.05  ~ "*",
#             TRUE             ~ "ns"
#           ),
#           label = paste0("p=", signif(p_value, 3), " ", sig)
#         )
#       
#       # ------------------------------
#       # 5. Annotation positions
#       # ------------------------------
#       y_max <- df %>%
#         dplyr::group_by(gene) %>%
#         dplyr::summarise(y = max(expression, na.rm = TRUE) * 1.15)
#       
#       annotation_df <- dplyr::left_join(pvals, y_max, by = "gene") %>%
#         dplyr::rename(p = label)
#       
#       # ------------------------------
#       # 6. Create boxplots for this dataset
#       # ------------------------------
#       p <- ggplot(df, aes(x = condition, y = expression, fill = condition)) +
#         geom_boxplot(outlier.shape = NA, alpha = 0.7) +
#         geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
#         facet_wrap(~ gene, scales = "free_y") +
#         geom_text(
#           data = annotation_df,
#           aes(x = 1.5, y = y, label = p),
#           inherit.aes = FALSE,
#           size = 3.5
#         ) +
#         labs(
#           title = paste0("Expression Boxplots: ", dataset_name),
#           x = NULL,
#           y = "Normalized Expression"  # changed 'Z-scored Expression' to 'Normalized Expression'
#         ) +
#         scale_fill_manual(values = c("NoNEC" = "#2CA02C", # chaning NoSepsis to NoNEC
#                                      "NEC"   = "#9467BD")) + # changing Sepsis to NEC
#         theme_bw(base_size = 12) +
#         theme(panel.grid = element_blank(),
#               legend.position = "none",
#               strip.background = element_blank(),
#               strip.text = element_text(face = "bold"))
#       
#       # ------------------------------
#       # 7. Optionally save plot
#       # ------------------------------
#       if (save_plots) {
#         fname <- paste0(prefix, "_", dataset_name, ".png")
#         ggsave(fname, plot = p, dpi = dpi, width = width, height = height, units = "in")
#       }
#       
#       # ------------------------------
#       # 8. Store results
#       # ------------------------------
#       results_list[[dataset_name]] <- list(
#         plot = p,
#         p_values = annotation_df,
#         data_used = df
#       )
#     }
#     
#     return(results_list)
#   }

############################################
# ---T-test--- 
############################################ 




############################################ 
# ---Differential Expression Analysis (limma)---
############################################

run_limma_DE <- function(expr_mat, pdata, condition_col = "condition",
                         lfc_threshold = NULL, sig_p = 0.05) {
  
  
  gsm_ids <- as.character(pdata$gsm)
  
  colnames(expr_mat) <- gsm_ids
  
  
  #-------------------------------
  #  Check that condition has ≥2 levels
  #-------------------------------
  cond_levels <- unique(as.character(pdata[[condition_col]]))
  if (length(cond_levels) < 2) {
    stop(paste0("Condition column '", condition_col, 
                "' has fewer than 2 levels: ", paste(cond_levels, collapse = ", ")))
  }
  
  #-------------------------------
  #  Build model matrix
  #-------------------------------
  design <- model.matrix(~ 0 + factor(pdata[[condition_col]]))
  
  # "0 +" removes the intercept so each level has its own column
  colnames(design) <- levels(factor(pdata[[condition_col]]))
  
  #-------------------------------
  #  Fit limma linear model and contrast
  #-------------------------------
  # Make sure condition is a factor
  cond <- factor(pdata[[condition_col]])
  if (length(levels(cond)) < 2) {
    stop(paste0("Condition column '", condition_col, "' has fewer than 2 levels"))
  }
  
  # Build design matrix: only the condition column, no intercept
  design <- model.matrix(~ 0 + cond)
  colnames(design) <- levels(cond)  # "NoSepsis" "Sepsis", etc.
  
  # Fit linear model
  fit <- limma::lmFit(expr_mat, design)
  
  # Create contrast: second level minus first level
  contrast_matrix <- limma::makeContrasts(
    contrasts = paste0(levels(cond)[2], "-", levels(cond)[1]),
    levels = design
  )
  
  # Apply contrast and compute statistics
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)
  
  # Extract full top table
  tt <- limma::topTable(fit2, number = nrow(expr_mat), adjust.method = "BH", sort.by = "P")
  
  if(nrow(tt) == 0){
    warning("topTable returned 0 rows. Returning empty table.")
    return(tt)
  }
  
  #-------------------------------
  #  Add SYMBOL column for mapping
  #-------------------------------
  tt$SYMBOL <- rownames(tt)
  
  #-------------------------------
  #  Manually add Expression column
  #-------------------------------
  if (is.null(lfc_threshold)) {
    q <- quantile(tt$logFC, na.rm = TRUE)
    up_cut <- q[4]     # 75%
    down_cut <- q[2]   # 25%
  } else {
    up_cut <- lfc_threshold
    down_cut <- -lfc_threshold
  }
  
  tt$Expression <- "Unchanged"
  tt$Expression[tt$logFC >= up_cut] <- "Up-regulated"
  tt$Expression[tt$logFC <= down_cut] <- "Down-regulated"
  
  #-------------------------------
  #  Manually add Significance column
  #-------------------------------
  tt$Significance <- "Non-significant"
  tt$Significance[tt$P.Value <= sig_p] <- "Significant"
  
  #-------------------------------
  #  Manually add Entrez IDs
  #-------------------------------
  tt$entrez <- mapIds(org.Hs.eg.db,
                      keys = tt$SYMBOL,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")
  
  #-------------------------------
  #  Optional: remove rows with NA entrez
  #-------------------------------
  tt <- tt[!is.na(tt$entrez), ]
  
  return(tt)
}


# Wrapper Function
run_limma_DE_list <- function(list_of_expr, list_of_pData, condition_col = "condition", genes_of_interest = NULL) {
  results_list <- list()
  skipped <- c()
  
  for (nm in names(list_of_expr)) {
    
    expr_mat <- list_of_expr[[nm]]
    
    # Filter matrix
    if (!is.null(genes_of_interest)) {
      expr_mat <- expr_mat[genes_of_interest, ]
    } else {
      expr_mat
    }
    
    pdata    <- list_of_pData[[nm]]
    
    # -------------------------------
    # 1️⃣ SAFELY align pdata with expression matrix
    # -------------------------------
    if ("gsm" %in% colnames(pdata)) {
      rownames(pdata) <- pdata$gsm
      pdata <- pdata[colnames(expr_mat), , drop = FALSE]   # safe reordering
      colnames(expr_mat) <- pdata$gsm
    } else {
      stop(paste("Dataset", nm, "has no 'gsm' column in pData"))
    }
    
    # -------------------------------
    # 2️⃣ CHECK for NA in condition column
    # -------------------------------
    if (any(is.na(pdata[[condition_col]]))) {
      message("Skipping dataset ", nm, " because condition column contains NA")
      skipped <- c(skipped, nm)
      next
    }
    
    # -------------------------------
    # 3️⃣ CHECK that condition has ≥2 levels
    # -------------------------------
    cond <- factor(pdata[[condition_col]], levels = unique(pdata[[condition_col]]))
    if (length(levels(cond)) < 2) {
      message("Skipping dataset ", nm,
              " because '", condition_col, "' has fewer than 2 levels: ",
              paste(levels(cond), collapse = ", "))
      skipped <- c(skipped, nm)
      next
    }
    
    # -------------------------------
    # 4️⃣ CALL run_limma_DE
    # -------------------------------
    res <- run_limma_DE(expr_mat = expr_mat,
                        pdata = pdata,
                        condition_col = condition_col)
    
    results_list[[nm]] <- res
  }
  
  attr(results_list, "skipped_datasets") <- skipped
  return(results_list)
}

############################################ 
# ---Differential Expression Analysis (limma)---
############################################


############################################ 
# ---Differential Expression Analysis (DESeq2)---
############################################

run_DESeq2_DE_nested <- function(nested_counts,
                                 list_of_pData,
                                 condition_col = "condition",
                                 lfc_threshold = NULL,
                                 sig_p = 0.05,
                                 genes_of_interest = NULL) {
  
  results_list <- list()
  skipped_datasets <- c()
  
  for (study_name in names(nested_counts)) {
    
    # -------------------------------
    # 1️⃣ Extract count matrix
    # -------------------------------
    count_mat <- nested_counts[[study_name]][[1]]  # always the first element
    
    # Filter for genes_of_interest
    if (!is.null(genes_of_interest)) {
      count_mat <- count_mat[genes_of_interest, ]
    } else {
      count_mat 
    }
    
    # -------------------------------
    # 2️⃣ Get corresponding pData from the separate list
    # -------------------------------
    if (!(study_name %in% names(list_of_pData))) {
      message("Skipping ", study_name, ": no pData found in list_of_pData")
      skipped_datasets <- c(skipped_datasets, study_name)
      next
    }
    
    pData_sub <- list_of_pData[[study_name]]
    pData_sub <- as.data.frame(pData_sub)
    pData_sub$gsm <- as.character(pData_sub$gsm)
    
    # -------------------------------
    # 3️⃣ Ensure GSMs match
    # -------------------------------
    if (!all(colnames(count_mat) %in% pData_sub$gsm)) {
      missing_gsm <- setdiff(colnames(count_mat), pData_sub$gsm)
      message("Skipping ", study_name, " because these GSM IDs are missing in pData: ",
              paste(missing_gsm, collapse = ", "))
      skipped_datasets <- c(skipped_datasets, study_name)
      next
    }
    
    # Reorder pData to match count matrix columns
    pData_sub <- pData_sub[match(colnames(count_mat), pData_sub$gsm), , drop = FALSE]
    
    # -------------------------------
    # 4️⃣ Check condition column
    # -------------------------------
    if (any(is.na(pData_sub[[condition_col]]))) {
      message("Skipping ", study_name, " because condition column contains NA")
      skipped_datasets <- c(skipped_datasets, study_name)
      next
    }
    
    cond <- factor(pData_sub[[condition_col]])
    if (length(levels(cond)) < 2) {
      message("Skipping ", study_name, " because condition has fewer than 2 levels")
      skipped_datasets <- c(skipped_datasets, study_name)
      next
    }
    
    # -------------------------------
    # 5️⃣ Run DESeq2
    # -------------------------------
    dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                                  colData   = pData_sub,
                                  design    = as.formula(paste0("~", condition_col)))
    dds <- DESeq(dds)
    
    res <- results(dds, contrast = c(condition_col, levels(cond)[2], levels(cond)[1]))
    res_df <- as.data.frame(res)
    
    # -------------------------------
    # 6️⃣ Add Expression and Significance columns
    # -------------------------------
    if (is.null(lfc_threshold)) {
      q <- quantile(res_df$log2FoldChange, na.rm = TRUE)
      up_cut <- q[4]     # 75%
      down_cut <- q[2]   # 25%
    } else {
      up_cut <- lfc_threshold
      down_cut <- -lfc_threshold
    }
    
    res_df <- res_df %>%
      dplyr::mutate(
        Expression = dplyr::case_when(
          log2FoldChange >= up_cut  ~ "Up-regulated",
          log2FoldChange <= -down_cut ~ "Down-regulated",
          TRUE ~ "Unchanged"
        ),
        Significance = dplyr::case_when(
          padj <= sig_p ~ "Significant",
          padj > sig_p  ~ "Non-significant"
        ),
        SYMBOL = rownames(res_df),
        entrez = AnnotationDbi::mapIds(org.Hs.eg.db,
                                       keys = rownames(res_df),
                                       column = "ENTREZID",
                                       keytype = "SYMBOL",
                                       multiVals = "first")
      ) %>%
      dplyr::arrange(padj)
    
    res_df <- na.omit(res_df)
    
    results_list[[study_name]] <- res_df
  }
  
  attr(results_list, "skipped_datasets") <- skipped_datasets
  return(results_list)
}

############################################ 
# ---Differential Expression Analysis (DESeq2)---
############################################


############################################ 
# ---Pathway Enrichment---
############################################

prepare_pathway_df <- function(up_list, down_list, n_show = 5) {
  
  # Function to process each pathway enrichment dataframe
  process_pathway <- function(df, direction, n_show) {
    if (is.null(df)) return(NULL)  # skip if NULL
    df %>%
      arrange(p.adjust) %>%
      slice_head(n = n_show) %>%
      mutate(direction = direction,
             minus_log10_p = -log10(p.adjust),
             term = Description)
  }
  
  # Process up-regulated pathways
  up_df <- imap(up_list, ~ process_pathway(.x, direction = 1, n_show = n_show)) %>%
    compact() %>%  # remove NULLs
    bind_rows(.id = "study")
  
  # Process down-regulated pathways
  down_df <- imap(down_list, ~ process_pathway(.x, direction = -1, n_show = n_show)) %>%
    compact() %>%  # remove NULLs
    bind_rows(.id = "study")
  
  # Combine up and down
  combined_df <- bind_rows(up_df, down_df) %>%
    mutate(count_dir = Count * direction,
           signed_logp = minus_log10_p * direction)
  
  return(combined_df)
}


# Plotting pathways
plot_pathways <- function(pathway_df, combine = TRUE, title = NULL, save_dir = NULL) {
  
  # Skip if pathway_df is empty
  if (nrow(pathway_df) == 0) {
    message("No pathways to plot. Exiting.")
    return(NULL)
  }
  
  # Determine save directory
  if (isTRUE(save_dir)) save_dir <- getwd()
  
  # Add signed logp and count_dir
  pathway_df <- pathway_df %>%
    dplyr::mutate(
      signed_logp = minus_log10_p * direction,
      count_dir = Count * direction
    ) %>%
    # Create term_ordered for plotting
    dplyr::arrange(direction, desc(abs(count_dir))) %>%
    dplyr::mutate(term_ordered = factor(term, levels = unique(term)))
  
  if (combine) {
    p <- ggplot(pathway_df, aes(x = term_ordered,
                                y = count_dir,
                                fill = signed_logp)) +
      geom_col(width = 0.6, colour = "black") +
      scale_fill_gradient2(low = "#7CAE00", mid = "white", high = "#C77CFF",
                           midpoint = 0, name = "-logP") +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
      coord_flip() +
      labs(x = NULL, y = "Gene counts", title = title) +
      theme_classic(base_size = 12) +
      theme(axis.text.y = element_text(hjust = 1),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            legend.title = element_text(size = 9),
            legend.position = "right",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
    
    # Save if directory is specified
    if (!is.null(save_dir)) {
      ggsave(filename = file.path(save_dir, "combined_pathways.png"),
             plot = p,
             width = 10,
             height = 6,
             dpi = 300,
             units = "in")
    }
    
    return(p)
    
  } else {
    # Split per study, skip studies with 0 rows
    plots <- pathway_df %>%
      split(.$study) %>%
      keep(~ nrow(.x) > 0) %>%  # skip empty studies
      map(~ ggplot(.x, aes(x = term_ordered,
                           y = count_dir,
                           fill = signed_logp)) +
            geom_col(width = 0.6, colour = "black") +
            scale_fill_gradient2(low = "#F68B1F", mid = "white", high = "#2166AC",
                                 midpoint = 0, name = "-logP") +
            geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
            coord_flip() +
            labs(x = NULL, y = "Gene counts", title = unique(.x$study)) +
            theme_classic(base_size = 12) +
            theme(axis.text.y = element_text(hjust = 1),
                  axis.ticks.y = element_blank(),
                  axis.line.y = element_blank(),
                  legend.title = element_text(size = 9),
                  legend.position = "right",
                  plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
      )
    
    # Save each plot if save_dir specified
    if (!is.null(save_dir)) {
      for (study_name in names(plots)) {
        ggsave(filename = file.path(save_dir, paste0(study_name, "_pathways.png")),
               plot = plots[[study_name]],
               width = 8,
               height = 6,
               dpi = 300,
               units = "in")
      }
    }
    
    return(plots)
  }
}

############################################ 
# ---Pathway Enrichment---
############################################


############################################ 
# ---Immune Cell Deconvolution---
############################################

plot_immune_deconv_list <- function(
    dna_list,
    rna_list,
    dna_pData,
    rna_pData,
    dna_method = "mcp_counter",
    rna_method = "quantiseq",
    cell_types = NULL,
    genes = NULL,                # <<---- NEW ARGUMENT
    save_path = NULL,
    dpi = 600
) {
  
  all_plots <- list()
  
  # ---- Function to process a single dataset ----
  process_dataset <- function(expr_matrix, pData, method) {
    
    # ---- OPTIONAL GENE FILTER ----
    if (!is.null(genes)) {
      keep <- intersect(genes, rownames(expr_matrix))
      if (length(keep) == 0) {
        stop("None of the provided genes were found in this expression matrix.")
      }
      expr_matrix <- expr_matrix[keep, , drop = FALSE]
    }
    # ----------------------------------------------
    
    # Deconvolution
    deconv_res <- immunedeconv::deconvolute(expr_matrix, method)
    
    # Remove rows with all NAs
    deconv_res <- deconv_res[rowSums(is.na(deconv_res[,-1])) < ncol(deconv_res[,-1]), ]
    
    # Pivot longer
    long_df <- deconv_res %>%
      tidyr::pivot_longer(
        cols = -cell_type,
        names_to = "gsm",
        values_to = "fraction"
      ) %>%
      dplyr::left_join(pData %>% dplyr::select(gsm, condition), by = "gsm")
    
    
    # ---- Normalize DNA MCP-counter scores to 0-1 ----
    if (method == "mcp_counter") {
      
      long_df <- long_df %>%
        dplyr::filter(!grepl("cytotoxic", cell_type, ignore.case = TRUE))
      
      long_df <- long_df %>%
        dplyr::group_by(cell_type) %>%
        dplyr::mutate(
          fraction = (fraction - min(fraction, na.rm = TRUE)) /
            (max(fraction, na.rm = TRUE) - min(fraction, na.rm = TRUE))
        ) %>%
        dplyr::ungroup()
    }
    
    # Summary
    long_df_summary <- long_df %>%
      dplyr::group_by(cell_type, condition) %>%
      dplyr::summarize(avg_fraction = mean(fraction), .groups = "drop")
    
    # P-values
    cell_pvals <- long_df %>%
      dplyr::group_by(cell_type) %>%
      dplyr::summarize(p_value = wilcox.test(fraction ~ condition)$p.value,
                       .groups = "drop") %>%
      dplyr::mutate(
        signif = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01  ~ "**",
          p_value < 0.05  ~ "*",
          TRUE            ~ "ns"
        ),
        label = paste0("p = ", signif(p_value, 3), " ", signif)
      )
    
    long_df_summary <- long_df_summary %>% left_join(cell_pvals, by = "cell_type")
    
    # Filter cell types
    if (!is.null(cell_types)) {
      filtered_cell_types <- intersect(cell_types, unique(long_df_summary$cell_type))
    } else {
      filtered_cell_types <- long_df_summary %>%
        dplyr::filter(avg_fraction > 0 & p_value <= 0.05) %>%
        dplyr::pull(cell_type) %>% unique()
    }
    
    long_df_summary <- long_df_summary %>% filter(cell_type %in% filtered_cell_types)
    
    pval_labels <- long_df_summary %>%
      dplyr::distinct(cell_type, label) %>%
      dplyr::mutate(y_pos = max(long_df_summary$avg_fraction) + 0.1)
    
    # Plot
    p <- ggplot2::ggplot(long_df_summary, aes(x = condition, y = avg_fraction, fill = condition)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), colour = "black", linewidth = 1.2) +
      ggplot2::geom_text(
        data = pval_labels,
        aes(x = 1.5, y = y_pos, label = label),
        inherit.aes = FALSE,
        size = 6
      ) +
      ggplot2::facet_wrap(~cell_type, scales = "free_y") +
      ggplot2::labs(
        x = "", y = "Average Cell Fraction", fill = "Condition"
      ) +
      ggplot2::scale_y_continuous(limits = c(0.0, 1)) +
      ggplot2::theme_classic() +
      ggplot2::scale_fill_manual(values = c("NoNEC" = "#2CA02C",
                                            "NEC"   = "#9467BD")) +
      ggplot2::theme(plot.margin = unit(c(.75,.75,.75,.75), "inches"),
                     strip.text = ggplot2::element_text(size = 14, face = "bold"),
                     axis.text = ggplot2::element_text(size = 16),
                     axis.title = ggplot2::element_text(size = 16),
                     legend.position = "none")
    
    return(list(plot = p, long_df = long_df, summary = long_df_summary, pvals = cell_pvals))
  }
  
  # ---- Updated workflow ----
  all_plots <- list()
  get_plot_name <- function(lst, prefix, i) {
    if (!is.null(names(lst)) && names(lst)[i] != "") names(lst)[i] else paste0(prefix, "_", i)
  }
  
  # DNA datasets
  for (i in seq_along(dna_list)) {
    expr_matrix <- dna_list[[i]]
    pData <- dna_pData[[i]]
    plot_name <- get_plot_name(dna_list, "DNA", i)
    all_plots[[plot_name]] <- process_dataset(expr_matrix, pData, dna_method)
  }
  
  # RNA datasets
  for (i in seq_along(rna_list)) {
    expr_matrix <- rna_list[[i]]
    pData <- rna_pData[[i]]
    plot_name <- get_plot_name(rna_list, "RNA", i)
    all_plots[[plot_name]] <- process_dataset(expr_matrix, pData, rna_method)
  }
  
  # Save images
  if (!is.null(save_path)) {
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    
    for (name in names(all_plots)) {
      ggsave(
        filename = file.path(save_path, paste0(name, ".png")),
        plot = all_plots[[name]]$plot,
        dpi = dpi, width = 10, height = 8, units = "in"
      )
    }
  }
  
  return(all_plots)
}


############################################ 
# ---Immune Cell Deconvolution---
############################################


############################################ 
# ---Removing invalid samples---
############################################ 
# expr_list <- matrices_dna[c("GSE25504", "GSE32472")]
# pdata_list <- sepsis_pData_dna[c("GSE25504", "GSE32472")]

clean_study_for_meta <- function(expr_list, pdata_list, id_col = "gsm") {
  cleaned_expr  <- list()
  cleaned_pdata <- list()
  
  for (study in names(expr_list)) {
    expr_mat <- expr_list[[study]]
    pdata_df <- pdata_list[[study]]
    
    if (!"condition" %in% colnames(pdata_df)) {
      stop("Study ", study, " does not have a 'condition' column")
    }
    if (!id_col %in% colnames(pdata_df)) {
      stop("Study ", study, " does not have an ID column '", id_col, "'")
    }
    
    # 1) Align by sample ID (GSM)
    sample_ids_expr <- colnames(expr_mat)
    sample_ids_pd   <- pdata_df[[id_col]]
    
    common_ids <- intersect(sample_ids_expr, sample_ids_pd)
    if (length(common_ids) == 0) {
      warning("Study ", study, " has no overlapping sample IDs between expr and pData. Skipping.")
      next
    }
    
    expr_mat <- expr_mat[, common_ids, drop = FALSE]
    pdata_df <- pdata_df[match(common_ids, sample_ids_pd), , drop = FALSE]
    
    if (ncol(expr_mat) != nrow(pdata_df)) {
      stop("After ID alignment in ", study, 
           " ncol(expr)=", ncol(expr_mat), 
           " != nrow(pData)=", nrow(pdata_df))
    }
    
    # 2) Clean condition (remove NAs / invalid)
    pdata_df$condition <- trimws(as.character(pdata_df$condition))
    
    valid_idx <- which(
      !is.na(pdata_df$condition) &
        pdata_df$condition != "" &
        pdata_df$condition != "#N/A!"
    )
    
    if (length(valid_idx) == 0) {
      warning("Study ", study, " has no valid condition labels after cleaning. Skipping.")
      next
    }
    
    pdata_df <- pdata_df[valid_idx, , drop = FALSE]
    expr_mat <- expr_mat[, valid_idx, drop = FALSE]
    
    if (ncol(expr_mat) != nrow(pdata_df)) {
      stop("After condition cleaning in ", study, 
           " ncol(expr)=", ncol(expr_mat), 
           " != nrow(pData)=", nrow(pdata_df))
    }
    
    if (length(unique(pdata_df$condition)) < 2) {
      warning("Study ", study, " does not contain both Case and Control after cleaning. Skipping.")
      next
    }
    
    pdata_df$condition <- factor(pdata_df$condition, levels = c("Control", "Case"))
    
    cleaned_expr[[study]]  <- expr_mat
    cleaned_pdata[[study]] <- pdata_df
  }
  
  list(
    expr_list  = cleaned_expr,
    pdata_list = cleaned_pdata
  )
}

############################################ 
# ---Removing invalid samples---
############################################



############################################ 
# ---Filtering for robust genes---
############################################

identify_robust_genes <- function(pool, fdr, p.het, g, lower_CI, upper_CI,
                                  effect_thresh, fdr_thresh, het_thresh) {
  robust_genes <- character()
  for (gene in rownames(g)) {
    g_gene <- g[gene, ]
    
    strong_effect <- abs(pool[gene]) > effect_thresh
    ci_consistent <- (lower_CI[gene] > 0) || (upper_CI[gene] < 0)
    hetero_ok     <- is.na(p.het[gene]) || p.het[gene] > het_thresh
    fdr_ok        <- fdr[gene] < fdr_thresh
    same_dir      <- (all(g_gene > 0)) || (all(g_gene < 0))
    
    if (strong_effect && ci_consistent && hetero_ok && fdr_ok && same_dir) {
      robust_genes <- c(robust_genes, gene)
    }
  }
  robust_genes
}

############################################ 
# ---Filtering for robust genes---
############################################


############################################ 
# ---Checking for overlapping GSMs---
############################################ 

check_gsm_overlap <- function(list_of_pData, gsm_col = "gsm") {
  
  gsm_df_list <- lapply(names(list_of_pData), function(study) {
    df <- list_of_pData[[study]]
    
    # Missing GSM column
    if (!gsm_col %in% colnames(df)) {
      warning(sprintf("Study %s has no '%s' column — skipped", study, gsm_col))
      return(NULL)
    }
    
    gsm <- as.character(df[[gsm_col]])
    
    # Empty GSM column
    if (length(gsm) == 0 || all(is.na(gsm))) {
      warning(sprintf("Study %s has empty GSMs — skipped", study))
      return(NULL)
    }
    
    data.frame(
      study = study,
      gsm   = gsm,
      stringsAsFactors = FALSE
    )
  })
  
  # Drop NULL entries
  gsm_df <- do.call(rbind, gsm_df_list[!sapply(gsm_df_list, is.null)])
  
  # Find overlaps
  overlaps <- unique(gsm_df$gsm[duplicated(gsm_df$gsm)])
  
  if (length(overlaps) == 0) {
    message("✅ No overlapping GSMs found across studies.")
    return(NULL)
  }
  
  gsm_df[gsm_df$gsm %in% overlaps, ]
}

############################################ 
# Checking for overlapping GSMs
############################################ 

# Usage
# check_gsm_overlap(sepsis_pData_dna)
# No overlapping GSMs


############################################
# ---Labeled PCA (LDA - Linear Discriminant Analysis)---
############################################
build_multicohort_dataset <- function(expr_list,
                                      pdata_list,
                                      genes,  # REQUIRED: robust genes from meta-analysis
                                      condition_col = "condition") {
  
  studies <- names(expr_list)
  names(pdata_list) <- studies
  
  
  # --------------------------
  # Subset expression matrices to robust genes
  # --------------------------
  expr_list <- lapply(expr_list, function(m){
    # Only keep genes in 'genes' that exist in this matrix
    m[intersect(genes, rownames(m)), , drop = FALSE]
  })
  
  # --------------------------
  # Combine matrices
  # --------------------------
  combined_matrix <- do.call(cbind, expr_list)   # genes x samples
  
  # --------------------------
  # Batch labels
  # --------------------------
  batch <- unlist(
    mapply(function(m, study){
      rep(study, ncol(m))
    }, expr_list, studies, SIMPLIFY = FALSE)
  )
  
  # --------------------------
  # Build metadata
  # --------------------------
  metadata <- do.call(rbind, lapply(seq_along(pdata_list), function(i){
    p <- pdata_list[[i]]
    data.frame(
      response = p[[condition_col]],
      study = studies[i],
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(metadata) <- colnames(combined_matrix)
  
  # --------------------------
  # Batch correction
  # --------------------------
  design <- model.matrix(~ response, data = metadata)
  batch_corrected <- limma::removeBatchEffect(
    combined_matrix,
    batch = batch,
    design = design
  )
  
  # --------------------------
  # ML-ready dataframe
  # --------------------------
  features <- t(batch_corrected)   # samples x genes
  multicohort_df <- cbind(metadata, features)
  
  # --------------------------
  # Return
  # --------------------------
  list(
    expression = batch_corrected,
    features = features,
    metadata = metadata,
    data = multicohort_df,
    batch = batch
  )
}


# features <- multicohort_df2 %>%
#   dplyr::select(-response, -study) %>% 
#   as.matrix()
# 
# labels <- case_when(multicohort_df2$response == "Control" ~ "Healthy Ctrl",
#                     multicohort_df2$response == "Case1" ~ "BPD",
#                     multicohort_df2$response == "Case2" ~ "Sepsis")
# 
# features_scaled <- scale(features)
# lda_res <- lda(labels ~ ., data = data.frame(features_scaled, labels = labels))
# 
# # Project the data onto LDA axes
# lda_pred <- predict(lda_res)$x
# lda_df <- data.frame(lda_pred, response = labels, study = multicohort_df2$study)


# Plotting
# png("\\\\ifs.win.uthscsa.edu/M1509-AhujaS/MainShare/Lois/Lois_Local/Dr_M/Projects/Neonatal_Sepsis/Figures/Prediction_Model/Multiclassification_ROC_curves/PCA_validation.png",
#     height = 8,
#     width = 8,
#     units = "in",
#     res = 800)
# ggplot(lda_df, aes(x = LD1, y = LD2, color = response)) +
#   geom_point(size = 3, alpha = 0.8) +
#   theme_classic(base_size = 14) +
#   theme(plot.margin = unit(c(1,1,1,1), "inches"),
#         panel.border = element_rect(color = "black", linewidth = 2),
#         legend.title = element_blank(),
#         legend.position = "inside",
#         legend.position.inside = c(.87, .1),
#         legend.background = element_rect(
#           fill = "white",
#           color = "black",
#           linewidth = 0.5,
#           linetype = "solid"
#         ),
#         legend.key = element_rect(fill = "white")) +
#   labs(x = "PC1",
#        y = "PC2",
#        color = "Condition") +
#   scale_color_manual(values = c("Healthy Ctrl" = "green3", "BPD" = "mediumpurple", "Sepsis" = "sienna3"))
# dev.off()

############################################
# ---Labeled PCA (LDA - Linear Discriminant Analysis)---
############################################
