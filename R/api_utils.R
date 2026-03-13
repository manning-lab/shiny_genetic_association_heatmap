# R/api_utils.R
# API utility functions for querying genetic association data from:
#   - GWAS Catalog REST API v2  (https://www.ebi.ac.uk/gwas/rest/api/v2)
#   - HuGeAMP BioIndex REST API (https://bioindex.hugeamp.org/api/bio)

# ---- Helpers -----------------------------------------------------------------

#' Null-coalescing operator: return first non-NULL / non-empty value
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b

#' Safely extract a character value from a possibly-nested list
safe_chr <- function(x, default = NA_character_) {
  tryCatch(as.character(x[[1]] %||% default), error = function(e) default)
}

#' Safely extract a numeric value
safe_num <- function(x, default = NA_real_) {
  tryCatch(suppressWarnings(as.numeric(x[[1]])) %||% default,
           error = function(e) default)
}

# ---- GWAS Catalog REST API v2 ------------------------------------------------

GWAS_BASE <- "https://www.ebi.ac.uk/gwas/rest/api/v2"

#' Fetch all associations for a given variant from GWAS Catalog API v2.
#'
#' @param variant_id rsID or variant ID (e.g. "rs7903146").
#' @param page_size  Maximum number of results per API page (default 500).
#' @return A data.frame with columns: variant_id, study_id, beta, pvalue,
#'         odds_ratio, trait, source.  Returns NULL on error / no data.
gwas_search_variant <- function(variant_id, page_size = 500) {
  url <- sprintf(
    "%s/associations?filter=variantId:%s&size=%d",
    GWAS_BASE,
    utils::URLencode(variant_id, reserved = FALSE),
    page_size
  )

  tryCatch({
    resp <- httr::GET(url,
                      httr::timeout(30),
                      httr::add_headers(Accept = "application/json"))

    if (httr::http_error(resp)) {
      message("GWAS API HTTP error ", httr::status_code(resp),
              " for variant: ", variant_id)
      return(NULL)
    }

    body <- httr::content(resp, "text", encoding = "UTF-8")
    data <- jsonlite::fromJSON(body, simplifyDataFrame = FALSE)

    assoc_list <- data[["_embedded"]][["associations"]]
    if (is.null(assoc_list) || length(assoc_list) == 0) return(NULL)

    rows <- lapply(assoc_list, function(a) {
      beta       <- safe_num(a[["beta"]])
      odds_ratio <- safe_num(a[["oddsRatio"]])
      pvalue     <- safe_num(a[["pValue"]])
      study_id   <- safe_chr(a[["study"]][["accessionId"]])
      var_id     <- safe_chr(a[["variant"]][["variantId"]], default = variant_id)
      trait <- tryCatch({
        et <- a[["efoTraits"]]
        if (length(et) > 0) as.character(et[[1]][["trait"]]) else NA_character_
      }, error = function(e) NA_character_)

      data.frame(
        variant_id = var_id,
        study_id   = study_id,
        beta       = beta,
        odds_ratio = odds_ratio,
        pvalue     = pvalue,
        trait      = trait,
        source     = "GWAS",
        stringsAsFactors = FALSE
      )
    })

    result <- do.call(rbind, rows)
    result <- result[!is.na(result$study_id), ]
    if (nrow(result) == 0) return(NULL)
    result

  }, error = function(e) {
    message("gwas_search_variant error: ", conditionMessage(e))
    NULL
  })
}

#' Fetch the beta value for a specific variant × study pair from GWAS Catalog v2.
#'
#' @param variant_id rsID or variant ID.
#' @param study_id   GWAS Catalog accession ID (e.g. "GCST001234").
#' @return Numeric beta value, or NA_real_ if unavailable.
gwas_get_beta <- function(variant_id, study_id) {
  url <- sprintf(
    "%s/associations?filter=variantId:%s,studyId:%s&size=10",
    GWAS_BASE,
    utils::URLencode(variant_id, reserved = FALSE),
    utils::URLencode(study_id,   reserved = FALSE)
  )

  tryCatch({
    resp <- httr::GET(url,
                      httr::timeout(30),
                      httr::add_headers(Accept = "application/json"))

    if (httr::http_error(resp)) return(NA_real_)

    body <- httr::content(resp, "text", encoding = "UTF-8")
    data <- jsonlite::fromJSON(body, simplifyDataFrame = FALSE)

    assoc_list <- data[["_embedded"]][["associations"]]
    if (is.null(assoc_list) || length(assoc_list) == 0) return(NA_real_)

    betas <- vapply(assoc_list,
                    function(a) safe_num(a[["beta"]]),
                    numeric(1))
    betas <- betas[!is.na(betas)]
    if (length(betas) == 0) return(NA_real_)
    betas[[1]]

  }, error = function(e) NA_real_)
}

# ---- HuGeAMP BioIndex REST API -----------------------------------------------

HUGEAMP_BASE <- "https://bioindex.hugeamp.org/api/bio"

#' Fetch all associations for a given variant from HuGeAMP BioIndex.
#'
#' @param variant_id rsID or variant ID (e.g. "rs7903146").
#' @param limit      Maximum records to retrieve (default 1000).
#' @return A data.frame with columns: variant_id, study_id, beta, pvalue,
#'         odds_ratio, trait, source.  Returns NULL on error / no data.
hugeamp_search_variant <- function(variant_id, limit = 1000) {
  url <- sprintf(
    "%s/query/associations?q=%s&limit=%d",
    HUGEAMP_BASE,
    utils::URLencode(variant_id, reserved = TRUE),
    limit
  )

  tryCatch({
    resp <- httr::GET(url,
                      httr::timeout(30),
                      httr::add_headers(Accept = "application/json"))

    if (httr::http_error(resp)) {
      message("HuGeAMP API HTTP error ", httr::status_code(resp),
              " for variant: ", variant_id)
      return(NULL)
    }

    body <- httr::content(resp, "text", encoding = "UTF-8")
    data <- jsonlite::fromJSON(body, simplifyDataFrame = TRUE)

    if (is.null(data[["data"]])) return(NULL)

    df <- if (is.data.frame(data[["data"]])) {
      data[["data"]]
    } else {
      tryCatch(as.data.frame(do.call(rbind, data[["data"]])),
               error = function(e) NULL)
    }

    if (is.null(df) || nrow(df) == 0) return(NULL)

    cols <- names(df)

    # Flexible column detection to handle API variations
    study_col  <- Find(function(c) c %in% cols,
                       c("phenotype", "dataset", "study", "ancestry", "trait"))
    beta_col   <- Find(function(c) c %in% cols,
                       c("beta", "effect", "effectSize", "log_or"))
    pval_col   <- Find(function(c) c %in% cols,
                       c("pValue", "pvalue", "p_value", "p", "pval"))
    var_col    <- Find(function(c) c %in% cols,
                       c("varId", "variant", "rsId", "variantId", "var"))
    trait_col  <- Find(function(c) c %in% cols,
                       c("phenotype", "trait", "dataset", "description"))
    or_col     <- Find(function(c) c %in% cols,
                       c("oddsRatio", "odds_ratio", "or"))

    result <- data.frame(
      variant_id = if (!is.null(var_col))
                     as.character(df[[var_col]])
                   else
                     rep(variant_id, nrow(df)),
      study_id   = if (!is.null(study_col))
                     as.character(df[[study_col]])
                   else
                     NA_character_,
      beta       = if (!is.null(beta_col))
                     suppressWarnings(as.numeric(df[[beta_col]]))
                   else
                     NA_real_,
      odds_ratio = if (!is.null(or_col))
                     suppressWarnings(as.numeric(df[[or_col]]))
                   else
                     NA_real_,
      pvalue     = if (!is.null(pval_col))
                     suppressWarnings(as.numeric(df[[pval_col]]))
                   else
                     NA_real_,
      trait      = if (!is.null(trait_col))
                     as.character(df[[trait_col]])
                   else
                     NA_character_,
      source     = "HuGeAMP",
      stringsAsFactors = FALSE
    )

    result <- result[!is.na(result$study_id) & nzchar(result$study_id), ]
    if (nrow(result) == 0) return(NULL)
    result

  }, error = function(e) {
    message("hugeamp_search_variant error: ", conditionMessage(e))
    NULL
  })
}

#' Fetch the beta value for a specific variant × study pair from HuGeAMP BioIndex.
#'
#' @param variant_id rsID or variant ID.
#' @param study_id   Phenotype / dataset identifier used in HuGeAMP.
#' @return Numeric beta value, or NA_real_ if unavailable.
hugeamp_get_beta <- function(variant_id, study_id) {
  result <- hugeamp_search_variant(variant_id)
  if (is.null(result) || nrow(result) == 0) return(NA_real_)

  matching <- result[result$study_id == study_id, ]
  if (nrow(matching) == 0) return(NA_real_)

  betas <- matching$beta[!is.na(matching$beta)]
  if (length(betas) == 0) return(NA_real_)
  betas[[1]]
}

# ---- Unified dispatcher -------------------------------------------------------

#' Search a variant using the specified data source.
#'
#' @param variant_id rsID or variant ID.
#' @param source     Either "gwas" or "hugeamp".
#' @return data.frame of associations or NULL.
search_variant <- function(variant_id, source = c("gwas", "hugeamp")) {
  source <- match.arg(source)
  switch(source,
    gwas    = gwas_search_variant(variant_id),
    hugeamp = hugeamp_search_variant(variant_id)
  )
}

#' Fetch beta for a variant × study pair from the specified source.
#'
#' @param variant_id rsID or variant ID.
#' @param study_id   Study / phenotype identifier.
#' @param source     Either "gwas" or "hugeamp".
#' @return Numeric beta, or NA_real_.
get_beta <- function(variant_id, study_id, source = c("gwas", "hugeamp")) {
  source <- match.arg(source)
  switch(source,
    gwas    = gwas_get_beta(variant_id, study_id),
    hugeamp = hugeamp_get_beta(variant_id, study_id)
  )
}
