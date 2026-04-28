#### Get normalized RNA-seq count matrix ####
  get_norm_RNA_counts <- function(expr_matrix, pData) {
    
    if (!"condition" %in% colnames(pData))
      stop("There must be a 'condition' column in pData")
    
    # Replace Inf/-Inf and NA with 0
    expr_matrix[!is.finite(expr_matrix)] <- 0
    
    dds <- DESeqDataSetFromMatrix(countData = round(expr_matrix), colData = pData, design = ~ condition)
    dds_1 <- DESeq(dds)
    
    norm_mtx <- counts(dds_1, normalized = TRUE)
    
    png("boxplot.png",
        units = "in",
        width = 8,
        height = 8,
        res = 600)
    layout(matrix(1:2, nrow=2))  # 2 rows, 1 column
    boxplot(expr_matrix, outline = FALSE, col = "skyblue", main = "Not Normalized")
    boxplot(norm_mtx, outline = FALSE, col = "skyblue", main = "Normalized")
    dev.off()
    
    return(norm_mtx)
  }
#### Get normalized RNA-seq count matrix ####
