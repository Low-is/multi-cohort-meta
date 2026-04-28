#### Generating expression matrices ####
generate_exprs_mtx <- function(DNA = NULL, RNA = NULL, dna_studies = list(), rna_studies = list()) {
  
  #### Handling duplicated genes ####
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
  
  #### For DNA microarray data sets ####
  if (!is.null(DNA) && DNA == TRUE) {
    if (!is.list(dna_studies))
      stop("dna_studies needs to be a list")
    
    if (!all(sapply(dna_studies, is.character)))
      stop("Elements in dna_studies must be character GEO IDs, e.g. list(study1 = 'GSEXXXX')")
    
    list_of_geo_data <- lapply(dna_studies, function(file) {
      tryCatch({
        gse_obj <- GEOquery::getGEO(GEO = file, AnnotGPL = TRUE)
        if (inherits(gse_obj, "ExpressionSet")) {
          list(gse_obj)
        } else {
          gse_obj
        }
      }, error = function(e) {
        warning(paste("Error retrieving", file, ":", e$message))
        NULL
      })
    })
    
    names(list_of_geo_data) <- names(dna_studies)
    
    out_dna <- vector("list", length(dna_studies))
    names(out_dna) <- names(dna_studies)
    
    for (s in seq_along(list_of_geo_data)) {
      es_list <- list_of_geo_data[[s]]
      if (is.null(es_list)) next
      
      per_platform_mtx <- vector("list", length(es_list))
      
      for (k in seq_along(es_list)) {
        geo <- es_list[[k]]
        if (is.null(geo)) next
        
        mtx <- Biobase::exprs(geo)
        fData <- Biobase::fData(geo)
        
        gene_cols <- c("Gene symbol", "Gene Symbol", "Symbol", "ILMN_Gene", "symbol", "GB_ACC", "ID")
        present <- intersect(gene_cols, colnames(fData))
        
        if (length(present) == 0L) {
          warning(sprintf("No gene annotation columns for study %s platform %d; using rownames.",
                          names(dna_studies)[s], k))
          fData$Genes <- rownames(fData)
        } else {
          
          # Priority-based selection
          if ("Gene symbol" %in% present) {
            selected_col <- "Gene symbol"
          } else if ("Gene Symbol" %in% present) {
            selected_col <- "Gene Symbol"
          } else if ("Symbol" %in% present) {
            selected_col <- "Symbol"
          } else if ("ILMN_Gene" %in% present) {
            selected_col <- "ILMN_Gene"
          } else if ("symbol" %in% present) {
            selected_col <- "symbol"
          } else if ("GB_ACC" %in% present) {
            selected_col <- "GB_ACC"
          } else if ("ID" %in% present) {
            selected_col <- "ID"
          } else {
            selected_col <- rownames(fData)
          }
          
          ## ================================
          ## NEW BRANCH FOR "ID"
          ## ================================
          if (selected_col == "ID") {
            probe_ids <- as.character(fData[[selected_col]])
            
            # Try mapping via hta20transcriptcluster.db
            symbols <- tryCatch(
              AnnotationDbi::mapIds(
                x = hta20transcriptcluster.db,
                keys = probe_ids,
                column = "SYMBOL",
                keytype = "PROBEID",
                multiVals = "first"
              ),
              error = function(e) {
                warning(sprintf(
                  "mapIds(PROBEID→SYMBOL) failed for study %s platform %d: %s. Using raw ID.",
                  names(dna_studies)[s], k, e$message
                ))
                rep(NA_character_, length(probe_ids))
              }
            )
            
            # if no annotation found, fall back to raw IDs
            if (all(is.na(symbols))) {
              warning(sprintf(
                "ID column in study %s platform %d did not produce SYMBOLs; using raw IDs.",
                names(dna_studies)[s], k
              ))
              fData$Genes <- probe_ids
            } else {
              fData$Genes <- symbols
            }
            
            ## ================================ 
            ## EXISTING GB_ACC BRANCH
            ## ================================
          } else if (selected_col == "GB_ACC") {
            
            probe_ids <- as.character(fData[[selected_col]])
            probe_ids_base <- sub("\\.\\d+$", "", probe_ids)
            
            symbols <- tryCatch(
              AnnotationDbi::mapIds(
                x = org.Hs.eg.db,
                keys = probe_ids_base,
                column = "SYMBOL",
                keytype = "ACCNUM",
                multiVals = "first"
              ),
              error = function(e) {
                warning(sprintf(
                  "mapIds ACCNUM failed for study %s platform %d: %s. Using raw GB_ACC.",
                  names(dna_studies)[s], k, e$message
                ))
                rep(NA_character_, length(probe_ids))
              }
            )
            
            if (all(is.na(symbols))) {
              warning(sprintf(
                "GB_ACC in study %s platform %d does not map as ACCNUM; using raw GB_ACC.",
                names(dna_studies)[s], k
              ))
              fData$Genes <- probe_ids
            } else {
              fData$Genes <- symbols
            }
            
            ## ================================
            ## ALL OTHER ANNOTATION COLUMNS
            ## ================================
          } else {
            fData$Genes <- as.character(fData[[selected_col]])
          }
        }
        
        genes_vector <- as.character(fData$Genes)
        
        df <- data.frame(Genes = genes_vector, mtx, check.names = FALSE)
        dt <- data.table::as.data.table(df)
        
        mtx2 <- handling_duplicates(dt)
        
        gene_names <- rownames(mtx2)
        exprs_mtx <- apply(mtx2, 2, as.numeric)
        log_mtx <- log(exprs_mtx + 1)
        log_mtx[is.nan(log_mtx)] <- 0
        
        norm_mtx <- limma::normalizeBetweenArrays(log_mtx, method = "quantile")
        rownames(norm_mtx) <- gene_names
        per_platform_mtx[[k]] <- norm_mtx
      }
      
      per_platform_mtx <- Filter(Negate(is.null), per_platform_mtx)
      if (length(per_platform_mtx) == 0L) next
      
      if (length(per_platform_mtx) == 1L) {
        out_dna[[s]] <- per_platform_mtx[[1]]
      } else {
        gene_lists <- lapply(per_platform_mtx, rownames)
        common_genes <- Reduce(intersect, gene_lists)
        
        if (length(common_genes) == 0L) {
          warning(sprintf(
            "No overlapping genes found across platforms for study %s; returning first platform matrix.",
            names(dna_studies)[s]
          ))
          out_dna[[s]] <- per_platform_mtx[[1]]
        } else {
          per_platform_common <- lapply(per_platform_mtx, function(m) {
            m_common <- m[common_genes, , drop = FALSE]
            m_common[common_genes, , drop = FALSE]
          })
          combined_mtx <- do.call(cbind, per_platform_common)
          rownames(combined_mtx) <- common_genes
          out_dna[[s]] <- combined_mtx
        }
      }
      
      pdf(paste0("boxplot_", names(dna_studies)[s], ".pdf"), width = 8, height = 8)
      boxplot(out_dna[[s]], outline = FALSE,
              col = "skyblue",
              main = paste0("Log transformation/Quantile\n", names(dna_studies)[s]))
      dev.off()
    }
    return(out_dna)
  }
  
  
  #### For RNA-Seq data sets #### 
  if (!is.null(RNA) && RNA == TRUE) {
    if (!is.list(rna_studies))
      stop("rna_studies needs to be a list of GEO accession IDs or GEO ExpressionSets")
    
    get_rnaseq_matrix <- function(geo_id, outdir = tempdir()) {
      if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
      
      # --- Attempt RAW tar download ---
      prefix <- substr(geo_id, 4, 6)
      folder <- paste0("GSE", prefix, "nnn")
      ftp_url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/", folder, "/", 
                        geo_id, "/suppl/", geo_id, "_RAW.tar")
      tar_path <- file.path(outdir, paste0(geo_id, "_RAW.tar"))
      
      message("Trying RAW tar for ", geo_id, " ...")
      status <- try(system(paste("wget -O", shQuote(tar_path), shQuote(ftp_url))), silent = TRUE)
      have_raw_tar <- file.exists(tar_path) && file.info(tar_path)$size > 0
      
      # --- Fallback to GEOquery ---
      if (inherits(status, "try-error") || !have_raw_tar) {
        message("RAW tar not available; fetching supp files via GEOquery for ", geo_id, " ...")
        supp_dir <- file.path(outdir, geo_id)
        GEOquery::getGEOSuppFiles(geo_id, baseDir = outdir, makeDirectory = TRUE, fetch_files = TRUE)
        extract_dir <- supp_dir
      } else {
        extract_dir <- file.path(outdir, paste0(geo_id, "_RAW"))
        if (!dir.exists(extract_dir)) dir.create(extract_dir, recursive = TRUE)
        system(paste("tar -xvf", shQuote(tar_path), "-C", shQuote(extract_dir)))
      }
      
      # --- Extract nested tar files if any ---
      tar_files <- list.files(extract_dir, pattern = "\\.tar$", full.names = TRUE)
      if (length(tar_files) > 0) {
        message("Extracting ", length(tar_files), " nested tar file(s) ...")
        for (tf in tar_files) system(paste("tar -xvf", shQuote(tf), "-C", shQuote(extract_dir)))
      }
      
      # --- Decompress all .txt.gz files ---
      gz_files <- list.files(extract_dir, pattern = "\\.txt\\.gz$", full.names = TRUE, recursive = TRUE)
      if (length(gz_files) > 0) {
        message("Decompressing ", length(gz_files), " .txt.gz files ...")
        lapply(gz_files, function(f) R.utils::gunzip(f, remove = FALSE, overwrite = TRUE))
      }
      
      # --- Gather all txt files ---
      txt_files <- list.files(extract_dir, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
      if (length(txt_files) == 0) stop("No RNA-seq .txt files found for ", geo_id)
      
      message("Detected ", length(txt_files), " processed RNA-seq .txt files for ", geo_id)
      
      # --- **FILE TYPE DETECTION BY FILENAME: count.txt = ENSG, others = FPKM/symbols** ---
      count_files <- txt_files[grepl("count\\.txt$", basename(txt_files))]
      fpkm_files <- txt_files[!grepl("count\\.txt$", basename(txt_files))]
      
      message("Count files (ENSG): ", length(count_files))
      message("FPKM/symbol files: ", length(fpkm_files))
      
      # --- Process COUNT FILES (ENSG → symbols) ---
      if (length(count_files) > 0) {
        dfs_count <- lapply(count_files, function(f) {
          dt <- data.table::fread(f)
          colnames(dt) <- c("gene", sub("\\.txt$", "", basename(f)))
          return(dt)
        })
        expr_dt_count <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), dfs_count)
        
        # ENSG → SYMBOL mapping
        message("Mapping ENSG from count files to symbols...")
        id_no_ver <- sub("\\..*$", "", expr_dt_count$gene)
        expr_dt_count <- expr_dt_count[!is.na(gene) & gene != "" & grepl("^ENSG", id_no_ver), ]
        id_no_ver <- sub("\\..*$", "", expr_dt_count$gene)
        
        gene_symbols <- AnnotationDbi::mapIds(
          org.Hs.eg.db, keys = id_no_ver, column = "SYMBOL", 
          keytype = "ENSEMBL", multiVals = "first"
        )
        expr_dt_count$gene_symbol <- gene_symbols
        expr_dt_count <- expr_dt_count[!is.na(gene_symbol), ]
        expr_dt_count <- expr_dt_count[, lapply(.SD, mean, na.rm = TRUE),
                                       by = gene_symbol,
                                       .SDcols = setdiff(names(expr_dt_count), c("gene", "gene_symbol"))]
        
      } else {
        expr_dt_count <- NULL
      }
      
      # --- Process FPKM FILES (keep symbols as-is) ---
      if (length(fpkm_files) > 0) {
        dfs_fpkm <- lapply(fpkm_files, function(f) {
          dt <- data.table::fread(f)
          if ("FPKM" %in% colnames(dt)) {
            gene_col <- intersect(c("gene_short_name", "tracking_id", "gene_id"), colnames(dt))[1]
            if (is.null(gene_col)) gene_col <- colnames(dt)[1]
            dt_sub <- dt[, .(Expr = mean(get("FPKM"), na.rm = TRUE)), by = get(gene_col)]
            colnames(dt_sub) <- c("gene", sub("\\.txt$", "", basename(f)))
          } else if (ncol(dt) == 2) {
            colnames(dt) <- c("gene", sub("\\.txt$", "", basename(f)))
            dt_sub <- dt
          } else {
            stop("Unknown FPKM format: ", f)
          }
          return(dt_sub)
        })
        expr_dt_fpkm <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), dfs_fpkm)
        expr_dt_fpkm <- expr_dt_fpkm[!is.na(gene) & gene != "", ]
        expr_dt_fpkm <- expr_dt_fpkm[, lapply(.SD, mean, na.rm = TRUE),
                                     by = gene,
                                     .SDcols = setdiff(names(expr_dt_fpkm), "gene")]
      } else {
        expr_dt_fpkm <- NULL
      }
      
      # --- COMBINE: prefer count-derived symbols, fallback to FPKM symbols ---
      if (!is.null(expr_dt_count)) {
        message("Using ENSG→symbol matrix from count files")
        expr_df <- as.data.frame(expr_dt_count)
        rownames(expr_df) <- expr_df$gene_symbol
        expr_df$gene_symbol <- NULL
      } else if (!is.null(expr_dt_fpkm)) {
        message("Using symbol matrix from FPKM files")
        expr_df <- as.data.frame(expr_dt_fpkm)
        rownames(expr_df) <- expr_df$gene
        expr_df$gene <- NULL
      } else {
        stop("No valid data found")
      }
      
      expr_matrix <- as.matrix(expr_df)
      
      # --- Phenotype data ---
      gse <- GEOquery::getGEO(geo_id, GSEMatrix = TRUE)[[1]]
      pdata <- pData(gse)
      gsm_ids <- colnames(expr_matrix)
      gsm_ids_short <- sub("^(GSM[0-9]+).*", "\\1", gsm_ids)
      colnames(expr_matrix) <- gsm_ids_short
      pdata <- pdata[gsm_ids_short, , drop = FALSE]
      
      message("Done: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " samples.")
      return(list(expr = expr_matrix, pData = pdata))
    }
    
    list_of_RNA_mtx <- lapply(rna_studies, function(rna) get_rnaseq_matrix(rna))
    names(list_of_RNA_mtx) <- rna_studies
    return(list_of_RNA_mtx)
  }
}
####
####
