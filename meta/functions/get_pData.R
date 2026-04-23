coalesce_cols <- function(df, sources, target) {
  present <- sources[sources %in% names(df)]
  
  if (length(present) == 0L) return(df)
  
  x <- df[[present[1]]]
  if (length(present) > 1L) {
    for (nm in present[-1]) {
      idx <- is.na(x) | x == ""
      x[idx] <- df[[nm]][idx]
    }
  }
  df[[target]] <- x
  # drop originals except target (if target is one of them)
  drop_these <- setdiff(present, target)
  df[drop_these] <- NULL
  df
}


get_pData <- function(studies = list()) {
  if (is.null(names(studies))) stop("Input must be a named list")
  
  out <- vector("list", length(studies))
  
  sqlite_file <- "C:\\Users\\RANDOLPHL\\Downloads\\GEOmetadb.sqlite"
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_file)
  
  for (i in seq_along(studies)) {
    study_id <- studies[[i]]
    study_name <- names(studies)[i]
    
    # ---------- submission year from GEOmetadb ----------
    placeholders <- paste(sprintf("'%s'", study_id), collapse = ", ")
    query_gsm <- sprintf(
      "SELECT DISTINCT
         gsm.gsm,
         COALESCE(substr(gse.submission_date, 1, 4), 'NA') AS year
       FROM gsm
       LEFT JOIN gse_gsm ON gsm.gsm = gse_gsm.gsm
       LEFT JOIN gse ON gse_gsm.gse = gse.gse
       WHERE gse.gse IN (%s)",
      placeholders
    )
    gsm_year_info <- DBI::dbGetQuery(con, query_gsm)
    
    # ---------- get the GSE soft / GSM list ----------
    gse_soft <- GEOquery::getGEO(study_id, AnnotGPL = TRUE, GSEMatrix = FALSE)
    gsm_list <- GEOquery::GSMList(gse_soft)
    
    platforms <- vapply(gsm_list, function(g) GEOquery::Meta(g)$platform_id, character(1))
    n_plat <- length(unique(platforms))
    
    if (n_plat > 1) {
      # ---------- Build a base_df that mimics pData(eset) for multi-platform ----------
      # Extract title, source_name, platform_id, and characteristics_ch1 (list) per GSM
      gsm_names <- names(gsm_list)
      titles <- vapply(gsm_list, function(g) { h <- GEOquery::Meta(g); if (!is.null(h$title)) h$title else NA_character_ }, character(1))
      sources <- vapply(gsm_list, function(g) { h <- GEOquery::Meta(g); if (!is.null(h$source_name_ch1)) h$source_name_ch1 else NA_character_ }, character(1))
      char_list <- lapply(gsm_list, function(g) {
        h <- GEOquery::Meta(g)
        ch <- h$characteristics_ch1
        if (is.null(ch) || length(ch) == 0L) return(NA_character_)
        return(as.character(ch))
      })
      
      # find maximum number of characteristics across GSMs
      max_chars <- max(vapply(char_list, function(x) if (all(is.na(x))) 0L else length(x), integer(1)))
      # create data.frame rows: pad each char vector to max_chars with NA
      char_mat <- t(sapply(char_list, function(x) {
        if (all(is.na(x))) {
          rep(NA_character_, max_chars)
        } else {
          length(x) <- max_chars
          x[is.na(x)] <- NA_character_
          as.character(x)
        }
      }, USE.NAMES = FALSE))
      
      # name characteristic columns exactly like R's Biobase does: characteristics_ch1, characteristics_ch1.1, ...
      if (max_chars == 0) {
        char_col_names <- character(0)
      } else if (max_chars == 1) {
        char_col_names <- "characteristics_ch1"
      } else {
        char_col_names <- c("characteristics_ch1", paste0("characteristics_ch1.", seq_len(max_chars - 1)))
      }
      
      # assemble base_df
      if (max_chars > 0) {
        base_df <- data.frame(
          gsm = gsm_names,
          title = titles,
          source_name_ch1 = sources,
          platform_id = platforms,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        # bind char columns
        base_df <- cbind(base_df, as.data.frame(char_mat, stringsAsFactors = FALSE, check.names = FALSE))
        names(base_df)[(ncol(base_df) - max_chars + 1):ncol(base_df)] <- char_col_names
      } else {
        base_df <- data.frame(
          gsm = gsm_names,
          title = titles,
          source_name_ch1 = sources,
          platform_id = platforms,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
      
      # ensure column names match expected names in your downstream code (characteristics_ch1 pattern)
      # (they already do above)
      
      # ---------- Now run the same pipeline you used in the n_plat == 1 branch ----------
      char_cols <- grep("^characteristics_ch1", names(base_df), value = TRUE)
      
      if (length(char_cols) > 0L) {
        # go long to parse key/value (exact same pipeline)
        pdata_long <- base_df %>%
          dplyr::select(gsm, platform_id, dplyr::all_of(char_cols)) %>%
          tidyr::pivot_longer(
            cols = dplyr::all_of(char_cols),
            names_to = "char_field",
            values_to = "raw"
          ) %>%
          dplyr::filter(!is.na(raw), raw != "") %>%
          tidyr::separate(
            raw, into = c("key", "value"),
            sep = ":", fill = "right", extra = "merge"
          ) %>%
          dplyr::mutate(
            key   = stringr::str_trim(key),
            value = stringr::str_trim(value)
          )
        
        # make characteristic names syntactically valid and go wide
        pdata_wide <- pdata_long %>%
          dplyr::mutate(key = make.names(key)) %>%
          tidyr::pivot_wider(
            id_cols    = c(gsm, platform_id),
            names_from = key,
            values_from = value
          )
        
        # add title / source_name (note base_df has source_name_ch1 column for multi-platform)
        pdata_long <- pdata_wide %>%
          dplyr::left_join(
            base_df %>%
              dplyr::transmute(
                gsm,
                title = if ("title" %in% names(base_df)) title else NA_character_,
                source_name = if ("source_name_ch1" %in% names(base_df)) source_name_ch1 else NA_character_
              ),
            by = "gsm"
          )
      } else {
        pdata_long <- base_df %>%
          dplyr::transmute(
            gsm,
            title       = if ("title" %in% names(base_df)) title else NA_character_,
            source_name = if ("source_name_ch1" %in% names(base_df)) source_name_ch1 else NA_character_,
            platform_id
          )
      }
      
    } else {
      # ---------- EXACTLY your original single-platform branch ----------
      eset_list <- GEOquery::getGEO(study_id, AnnotGPL = TRUE, GSEMatrix = TRUE)
      eset <- eset_list[[1]]
      base_df <- as.data.frame(Biobase::pData(eset), stringsAsFactors = FALSE)
      base_df$gsm <- rownames(base_df)
      base_df$platform_id <- unique(platforms)
      
      char_cols <- grep("^characteristics_ch1", names(base_df), value = TRUE)
      
      if (length(char_cols) > 0L) {
        # go long to parse key/value
        pdata_long <- base_df %>%
          dplyr::select(gsm, platform_id, dplyr::all_of(char_cols)) %>%
          tidyr::pivot_longer(
            cols = dplyr::all_of(char_cols),
            names_to = "char_field",
            values_to = "raw"
          ) %>%
          dplyr::filter(!is.na(raw), raw != "") %>%
          tidyr::separate(
            raw, into = c("key", "value"),
            sep = ":", fill = "right", extra = "merge"
          ) %>%
          dplyr::mutate(
            key   = stringr::str_trim(key),
            value = stringr::str_trim(value)
          )
        
        # make characteristic names syntactically valid and go wide
        pdata_wide <- pdata_long %>%
          dplyr::mutate(key = make.names(key)) %>%
          tidyr::pivot_wider(
            id_cols    = c(gsm, platform_id),
            names_from = key,
            values_from = value
          )
        
        # add title / source_name
        pdata_long <- pdata_wide %>%
          dplyr::left_join(
            base_df %>%
              dplyr::transmute(
                gsm,
                title       = if ("title" %in% names(base_df)) title else NA_character_,
                source_name = if ("source_name_ch1" %in% names(base_df)) source_name_ch1 else NA_character_
              ),
            by = "gsm"
          )
      } else {
        pdata_long <- base_df %>%
          dplyr::transmute(
            gsm,
            title       = if ("title" %in% names(base_df)) title else NA_character_,
            source_name = if ("source_name_ch1" %in% names(base_df)) source_name_ch1 else NA_character_,
            platform_id
          )
      }
    } # end n_plat if/else
    
    # ---------- Common post-processing (same as you had) ----------
    pdata_long <- dplyr::left_join(pdata_long, gsm_year_info, by = "gsm")
    
    ga_sources <- c("corrected.gestational.age..weeks.days.",
                    "corrected.gestational.age.at.sample..weeks.days.")
    pdata_long <- coalesce_cols(pdata_long, sources = ga_sources, target = "corrected_gestational_age")
    
    pdata_long$study <- study_name
    out[[i]] <- pdata_long
  } # end studies loop
  
  names(out) <- names(studies)
  DBI::dbDisconnect(con)
  return(out)
}


add_condition_column <- function(df, column, case_patterns, control_patterns){ # This may be more challenging to implement
  
  if(!(column %in% colnames(df))) stop("Column not found in df")
  
  case_regex <- paste(case_patterns, collapse = "|")
  control_regex <- paste(control_patterns, collapse = "|")
  
  df$condition <- dplyr::case_when(
  grepl(case_regex, tolower(df[[column]])) ~ "Case",
  grepl(control_regex, tolower(df[[column]])) ~ "Control",
  TRUE ~ NA_character_
)
  
  df$condition <- factor(df$condition, levels = c("Control", "Case"), exclude = NULL)
  
  return(df)
}



detect_condition_column <- function(df, case_patterns, control_patterns) {

  normalize <- function(x) {
    x <- tolower(as.character(x))
    x <- gsub("[^a-z0-9 ]", " ", x)
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }

  char_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
  char_cols <- setdiff(char_cols, c("gsm", "study", "platform_id"))

  case_regex <- paste(normalize(case_patterns), collapse = "|")
  control_regex <- paste(normalize(control_patterns), collapse = "|")

  scores <- sapply(char_cols, function(col) {

    vals <- normalize(df[[col]])
    vals <- vals[!is.na(vals) & vals != ""]
    vals <- unique(vals)

    n_unique <- length(vals)
    score <- 0

    if (n_unique == 2) score <- score + 5
    if (n_unique <= 5) score <- score + 2

    case_hit <- any(grepl(case_regex, vals))
    control_hit <- any(grepl(control_regex, vals))

    if (case_hit && control_hit) {
      score <- score + 6
    } else if (case_hit || control_hit) {
      score <- score + 3
    }

    return(score)
  })

  if (length(scores) == 0) return(NULL)

  best <- names(which.max(scores))

  if (scores[best] == 0) return(NULL)

  return(best)
}



apply_condition_to_list <- function(pdata_list, case_patterns, control_patterns) {
  
  normalize <- function(x) {
    x <- tolower(as.character(x))
    x <- gsub("[^a-z0-9 ]", " ", x)
    x <- gsub("\\s+", " ", x)
    trimws(x)
  }
  
  out <- lapply(names(pdata_list), function(study) {
    
    df <- pdata_list[[study]]
    
    message("Processing: ", study)
    
    col <- detect_condition_column(
      df,
      case_patterns = case_patterns,
      control_patterns = control_patterns
    )
    
    if (is.null(col)) {
      warning("No condition column detected for ", study)
      return(NULL)
    }
    
    message("Using column: ", col)
    
    values <- normalize(df[[col]])
    
    case_regex <- paste(normalize(case_patterns), collapse = "|")
    control_regex <- paste(normalize(control_patterns), collapse = "|")
    
    df$condition <- dplyr::case_when(
      grepl(case_regex, values) ~ "Case",
      grepl(control_regex, values) ~ "Control",
      TRUE ~ NA_character_
    )
    
    if (all(is.na(df$condition))) {
      warning("Skipping study (no Case/Control signal): ", study)
      return(NULL)
    }
    
    df$condition <- factor(df$condition, levels = c("Control", "Case"))
    
    return(df)
  })
  
  # KEEP ONLY VALID ONES + FIX NAMES SAFELY
  names(out) <- names(pdata_list)
  out <- out[!sapply(out, is.null)]
  
  keep <- !sapply(out, is.null)
  out <- out[keep]
  
  names(out) <- names(pdata_list)[keep]
  
  return(out)
}
