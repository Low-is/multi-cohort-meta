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



resolve_matrix_names <- function(expr_colnames, expr_matrix, pdata) {

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

  # safety check
  if (length(expr_colnames) == 0) {
    warning("No valid sample columns after filtering")
    return(NULL)
  }

  # -------------------------
  # DETECT GSM MODE
  # -------------------------

  is_gsm <- mean(grepl("^GSM", expr_colnames)) > 0.8

  # -------------------------
  # CASE 1: GSM-BASED MATCHING
  # -------------------------

  if (is_gsm) {

    matched_rows <- match(expr_colnames, pdata$gsm)

    if (any(is.na(matched_rows))) {
      warning("GSM matching failed")
      return(NULL)
    }

    if (length(matched_rows) != length(expr_colnames)) {
      warning("GSM mismatch: length issue")
      return(NULL)
    }

    if (any(duplicated(matched_rows))) {
      warning("Duplicate GSM mapping detected")
      return(NULL)
    }

    message("Matched using GSM IDs")

    ord <- order(matched_rows)

    expr_matrix <- expr_matrix[, ord, drop = FALSE]
    colnames(expr_matrix) <- pdata$gsm[matched_rows][ord]

    return(expr_matrix)
  }

# -------------------------
  # CASE 2: TITLE-BASED MATCHING
  # -------------------------

  candidates <- generate_candidate_titles(as.character(pdata$title))

  expr_norm <- normalize(expr_colnames)

  for (method in names(candidates)) {

    candidate_norm <- normalize(candidates[[method]])

    # skip ambiguous mappings
    if (any(duplicated(candidate_norm))) {
      next
    }

    matched_rows <- match(expr_norm, candidate_norm)

    # keep only matched samples
    keep <- !is.na(matched_rows)

    if (sum(keep) == 0) {
      next
    }

    matched_rows <- matched_rows[keep]

    # subset matrix to matched samples only
    expr_subset <- expr_matrix[, keep, drop = FALSE]

    # avoid duplicate mappings
    if (any(duplicated(matched_rows))) {
      next
    }

    # reorder columns to pdata order
    ord <- order(matched_rows)

    matched_rows <- matched_rows[ord]

    expr_subset <- expr_subset[, ord, drop = FALSE]

    # rename columns
    colnames(expr_subset) <- pdata$gsm[matched_rows]

    message(
      "Matched using: ",
      method,
      " (",
      ncol(expr_subset),
      " samples)"
    )

    return(expr_subset)
  }

  warning("Could not resolve matrix sample names")
  return(NULL)

}
