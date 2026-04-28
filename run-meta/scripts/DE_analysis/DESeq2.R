# Perform meta-analysis (DESeq2)
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
    
    # Replace Inf/-Inf and NA with 0
    count_mat[!is.finite(count_mat)] <- 0
    
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
