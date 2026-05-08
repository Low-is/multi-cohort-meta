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

  # -------------------------
  # GENERATE TITLE CANDIDATES
  # -------------------------

  candidates <- generate_candidate_titles(
    as.character(pdata$title)
  )

  expr_norm <- normalize(expr_colnames)

  # -------------------------
  # TRY EACH MATCHING STRATEGY
  # -------------------------

  for (method in names(candidates)) {

    candidate <- candidates[[method]]

    candidate_norm <- normalize(candidate)

    # -------------------------
    # REQUIRE ALL MATRIX COLS TO MATCH
    # -------------------------

    if (all(expr_norm %in% candidate_norm)) {

      matched_rows <- match(
        expr_norm,
        candidate_norm
      )

      # -------------------------
      # VALIDATE MATCHES
      # -------------------------

      if (any(is.na(matched_rows))) {
        next
      }

      # avoid duplicated mappings
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

        matched_titles =
          candidate[matched_rows],

        matched_gsms =
          pdata$gsm[matched_rows]
      ))
    }
  }

  warning("Could not resolve matrix sample names")

  return(NULL)
}
