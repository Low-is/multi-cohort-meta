 handling_duplicates <- function(data_table) {
    if (is.list(data_table$Genes)) {
      data_table[, Genes := vapply(Genes, function(x) {
        if (is.null(x) || length(x) == 0 || all(is.na(x))) {
          return(NA_character_)
        } else if (is.list(x)) {
          return(as.character(unlist(x)))
        } else {
          return(as.character(x))
        }
      }, FUN.VALUE = character(1))]
    }
    
    data_table <- na.omit(data_table[(Genes != "" & !is.na(Genes)), ])
    
    if (any(duplicated(data_table$Genes))) {
      exprs_data.table <- data_table[, lapply(.SD, mean), by = Genes, .SDcols = !'Genes']
    } else {
      exprs_data.table <- data_table
    }
    
    exprs_mtx <- as.matrix(exprs_data.table[, -1])
    rownames(exprs_mtx) <- exprs_data.table$Genes
    
    return(exprs_mtx)
  }
