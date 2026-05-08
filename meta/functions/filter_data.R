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



resolve_matrix_names <- function(expr_colnames, pdata) {

  # -------------------------
  # REMOVE OBVIOUS NON-SAMPLES
  # -------------------------

  expr_colnames <- expr_colnames[
    !tolower(expr_colnames) %in% c(
      "id",
      "gene",
      "genes",
      "length",
      "symbol"
    )
  ]

  expr_norm <- normalize(expr_colnames)

  # -------------------------
  # DETECT GSM MODE
  # -------------------------

  is_gsm <- all(grepl("^GSM", expr_colnames))

  # -------------------------
  # CASE 1: GSM-BASED MATCHING
  # -------------------------

  if (is_gsm) {

    gsm_norm <- normalize(pdata$gsm)

    matched_rows <- match(expr_norm, gsm_norm)

    if (any(is.na(matched_rows))) {
      warning("GSM matching failed")
      return(NULL)
    }

    if (any(duplicated(matched_rows))) {
      warning("Duplicate GSM mapping detected")
      return(NULL)
    }

    message("Matched using GSM IDs")

    return(list(
      method = "gsm",
      pdata_rows = matched_rows,
      matched_titles = pdata$title[matched_rows],
      matched_gsms = pdata$gsm[matched_rows]
    ))
  }

  # -------------------------
  # CASE 2: TITLE-BASED MATCHING
  # -------------------------

  candidates <- generate_candidate_titles(
    as.character(pdata$title)
  )

  for (method in names(candidates)) {

    candidate <- candidates[[method]]
    candidate_norm <- normalize(candidate)

    if (all(expr_norm %in% candidate_norm)) {

      matched_rows <- match(expr_norm, candidate_norm)

      # -------------------------
      # VALIDATE MATCHES
      # -------------------------

      if (any(is.na(matched_rows))) {
        next
      }

      if (any(duplicated(matched_rows))) {
        warning(
          "Duplicate sample mapping detected using: ",
          method
        )
        next
      }

      message("Matched using: ", method)

      return(list(
        method = method,
        pdata_rows = matched_rows,
        matched_titles = candidate[matched_rows],
        matched_gsms = pdata$gsm[matched_rows]
      ))
    }
  }

  warning("Could not resolve matrix sample names")

  return(NULL)
}
