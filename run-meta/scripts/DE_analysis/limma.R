# Perform meta-analysis (Differential Analysis)

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
  tt <- tt[tt$SYMBOL %in% keys(org.Hs.eg.db, keytype = "SYMBOL"), ]
  
  if(nrow(tt) == 0){
    warning("No valid SYMBOLs left after filtering. Returning empty table.")
    return(tt)
  }
  
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
