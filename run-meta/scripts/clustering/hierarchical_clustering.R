plot_hclust_metagenes <- function(expr_list,
                                  pData_list,
                                  meta_genes,
                                  gsm_col = "gsm",        # column to use for rownames if needed
                                  save_path = NULL,
                                  prefix = "hclust_meta",
                                  dpi = 600,
                                  width = 8,
                                  height = 6) {
  # Load libraries
  require(dplyr)
  
  # Check list lengths
  if (length(expr_list) != length(pData_list)) {
    stop("expr_list and pData_list must have the same number of elements.")
  }
  
  results <- list()
  
  for (i in seq_along(expr_list)) {
    expr <- expr_list[[i]]
    pdata <- pData_list[[i]]
    
    # Dataset name
    dataset_name <- names(expr_list)[i]
    if (is.null(dataset_name) || dataset_name == "") {
      dataset_name <- paste0("dataset_", i)
    }
    
    # Ensure meta genes exist
    valid_genes <- intersect(meta_genes, rownames(expr))
    if (length(valid_genes) == 0) {
      warning(paste("No meta genes found in", dataset_name, "- skipping"))
      next
    }
    
    # Subset expression matrix
    sub_expr <- expr[valid_genes, , drop = FALSE]
    
    # Ensure condition column exists
    if (!"condition" %in% colnames(pdata)) {
      stop(paste("pData for", dataset_name, "must contain column 'condition'"))
    }
    
    # Fix rownames of pData if necessary
    if (!all(colnames(sub_expr) %in% rownames(pdata))) {
      if (gsm_col %in% colnames(pdata)) {
        rownames(pdata) <- as.character(pdata[[gsm_col]])
      } else {
        stop(paste("pData rownames do not match expression matrix columns and column", gsm_col, "not found for", dataset_name))
      }
    }
    
    # Reorder pData to match columns
    pdata <- pdata[colnames(sub_expr), , drop = FALSE]
    
    # Z-score scaling of genes
    scaled_expr <- t(scale(t(sub_expr)))  # genes × samples
    
    # Distance and hierarchical clustering
    dist_mat <- dist(t(scaled_expr), method = "euclidean")  # distance between samples
    hc <- hclust(dist_mat, method = "ward.D2")
    
    # Prepare file path
    out_file <- paste0(prefix, "_", dataset_name, ".png")
    if (!is.null(save_path)) {
      if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
      out_file <- file.path(save_path, out_file)
    }
    
    # Save dendrogram as PNG
    png(filename = out_file, width = width, height = height, units = "in", res = dpi)
    plot(hc, main = paste0("Hierarchical Clustering: ", dataset_name),
         xlab = "", sub = "", cex = 0.8)
    dev.off()
    
    results[[dataset_name]] <- list(
      hc = hc,
      scaled_matrix = scaled_expr,
      pdata = pdata,
      file = out_file
    )
  }
  
  return(results)
}
