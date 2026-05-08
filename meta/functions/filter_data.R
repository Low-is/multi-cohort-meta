generate_candidate_titles <- function(title_vec) {

  list(

    # original
    full = title_vec,

    # remove after comma
    before_comma =
      sub(",.*$", "", title_vec),

    # remove after semicolon
    before_semicolon =
      sub(";.*$", "", title_vec),

    # remove after whitespace
    first_token =
      sub(" .*$", "", title_vec),

    # remove trailing underscore suffix
    before_underscore_suffix =
      sub("_[^_]+$", "", title_vec),

    # remove trailing dash suffix
    before_dash_suffix =
      sub("-[^-]+$", "", title_vec),

    # normalized versions
    normalized =
      normalize(title_vec),

    normalized_before_comma =
      normalize(sub(",.*$", "", title_vec)),

    normalized_before_underscore =
      normalize(sub("_[^_]+$", "", title_vec)),

    normalized_before_dash =
      normalize(sub("-[^-]+$", "", title_vec))
  )
}



resolve_matrix_names(expr_colnames, expr_matrix, pdata) {

  # -------------------------
  # REMOVE OBVIOUS NON-SAMPLES
  # -------------------------

  keep <- !tolower(expr_colnames) %in% c(
    "id",
    "gene",
    "genes",
    "length",
    "symbol"
  )

  expr_colnames <- expr_colnames[keep]
  expr_matrix <- expr_matrix[, keep, drop = FALSE]

  # -------------------------
  # DETECT GSM MODE
  # -------------------------

  is_gsm <- all(grepl("^GSM", expr_colnames))

  # -------------------------
  # CASE 1: GSM-BASED MATCHING
  # -------------------------

  if (is_gsm) {

    matched_rows <- match(expr_colnames, pdata$gsm)

    if (any(is.na(matched_rows))) {
      warning("GSM matching failed")
      return(NULL)
    }

    if (any(duplicated(matched_rows))) {
      warning("Duplicate GSM mapping detected")
      return(NULL)
    }

    message("Matched using GSM IDs")

    expr_matrix <- expr_matrix[, order(matched_rows), drop = FALSE]
    colnames(expr_matrix) <- pdata$gsm[matched_rows][order(matched_rows)]

    return(expr_matrix)
  }
  
  # -------------------------
  # CASE 2: TITLE-BASED MATCHING
  # -------------------------

  candidates <- generate_candidate_titles(as.character(pdata$title))

  for (method in names(candidates)) {

    candidate <- candidates[[method]]
    candidate_norm <- normalize(candidate)

    expr_norm <- normalize(expr_colnames)

    if (all(expr_norm %in% candidate_norm)) {

      matched_rows <- match(expr_norm, candidate_norm)

      if (any(is.na(matched_rows))) next
      if (any(duplicated(matched_rows))) next

      message("Matched using: ", method)

      expr_matrix <- expr_matrix[, order(matched_rows), drop = FALSE]
      colnames(expr_matrix) <- pdata$gsm[matched_rows][order(matched_rows)]

      return(expr_matrix)
    }
  }

  warning("Could not resolve matrix sample names")
  return(NULL)
}
