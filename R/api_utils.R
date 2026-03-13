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

# Pattern for recognising rsID-style variant identifiers
RSID_PATTERN <- "^rs[0-9]+$"

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

# Ensembl REST API base used for rsID → GRCh38 position resolution
ENSEMBL_BASE <- "https://rest.ensembl.org"

#' Resolve an rsID to a GRCh38 chromosomal region using the Ensembl REST API.
#'
#' @param rsid  An rsID string, e.g. "rs979012".
#' @param flank Bases added on each side of the SNP position for a region query
#'   (default 0 — returns the exact SNP region; increase if you need a window).
#' @return A list with elements \code{chrom} (character), \code{start} (integer),
#'   and \code{end} (integer) on GRCh38, or NULL if resolution fails.
resolve_rsid_to_grch38 <- function(rsid, flank = 0L) {
  url <- sprintf(
    "%s/variation/human/%s?content-type=application%%2Fjson",
    ENSEMBL_BASE,
    utils::URLencode(rsid, reserved = FALSE)
  )

  tryCatch({
    resp <- httr::GET(url,
                      httr::timeout(20),
                      httr::add_headers(Accept = "application/json"))

    if (httr::http_error(resp)) {
      message("Ensembl API HTTP error ", httr::status_code(resp),
              " for rsID: ", rsid)
      return(NULL)
    }

    body <- httr::content(resp, "text", encoding = "UTF-8")
    data <- jsonlite::fromJSON(body, simplifyDataFrame = FALSE)

    mappings <- data[["mappings"]]
    if (is.null(mappings) || length(mappings) == 0) return(NULL)

    # Prefer GRCh38 assembly, otherwise take the first mapping
    grch38 <- Filter(
      function(m) identical(m[["assembly_name"]], "GRCh38"),
      mappings
    )
    m <- if (length(grch38) > 0) grch38[[1]] else mappings[[1]]

    chrom <- as.character(m[["seq_region_name"]])
    start <- as.integer(m[["start"]])
    end   <- as.integer(m[["end"]])

    if (is.na(chrom) || is.na(start)) return(NULL)

    list(chrom = chrom,
         start = max(1L, start - flank),
         end   = end + flank)

  }, error = function(e) {
    message("resolve_rsid_to_grch38 error: ", conditionMessage(e))
    NULL
  })
}

#' Parse a raw HuGeAMP BioIndex API response body into a tidy data.frame.
#'
#' Used internally by both the rsID and position-based query paths.
#'
#' @param body         Raw response body string.
#' @param variant_id   Fallback variant identifier for the \code{variant_id}
#'   column when the API response does not contain one.
#' @return A data.frame or NULL.
.parse_hugeamp_body <- function(body, variant_id) {
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
}

#' Fetch all associations for a given variant from HuGeAMP BioIndex.
#'
#' Queries the BioIndex first by rsID.  If the server returns HTTP 400
#' (malformed query — common when the variant is not indexed by rsID), the
#' function automatically resolves the rsID to a GRCh38 chromosomal region via
#' the Ensembl REST API and retries using a positional query.
#'
#' @param variant_id rsID or variant ID (e.g. "rs7903146").
#' @param limit      Maximum records to retrieve (default 1000).
#' @return A data.frame with columns: variant_id, study_id, beta, pvalue,
#'         odds_ratio, trait, source.  Returns NULL on error / no data.
hugeamp_search_variant <- function(variant_id, limit = 1000) {
  url <- build_hugeamp_query_url(variant_id, limit = limit)

  tryCatch({
    resp <- httr::GET(url,
                      httr::timeout(30),
                      httr::add_headers(Accept = "application/json"))

    status <- httr::status_code(resp)

    # ---- Happy path ----------------------------------------------------------
    if (!httr::http_error(resp)) {
      body <- httr::content(resp, "text", encoding = "UTF-8")
      return(.parse_hugeamp_body(body, variant_id))
    }

    # ---- Fallback: resolve rsID to GRCh38 position and retry ----------------
    # A 400 typically means the BioIndex does not recognise the rsID as a valid
    # query key.  Many variants are indexed by chromosomal varId
    # (chr:pos:ref:alt) rather than rsID; in that case we look up the position
    # via Ensembl and query by region.
    if (status == 400L && grepl(RSID_PATTERN, variant_id, ignore.case = TRUE)) {
      message("HuGeAMP API returned 400 for rsID '", variant_id,
              "'. Trying positional fallback via Ensembl …")

      coords <- resolve_rsid_to_grch38(variant_id)

      if (!is.null(coords)) {
        pos_url <- build_hugeamp_query_url(
          sprintf("%s:%d-%d", coords$chrom, coords$start, coords$end),
          limit = limit
        )
        resp2 <- httr::GET(pos_url,
                           httr::timeout(30),
                           httr::add_headers(Accept = "application/json"))

        if (!httr::http_error(resp2)) {
          body2 <- httr::content(resp2, "text", encoding = "UTF-8")
          result <- .parse_hugeamp_body(body2, variant_id)
          if (!is.null(result) && nrow(result) > 0) {
            message("Positional fallback succeeded for '", variant_id, "'.")
            return(result)
          }
        }
      }
    }

    message("HuGeAMP API HTTP error ", status, " for variant: ", variant_id)
    NULL

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

# ---- Query URL builders -------------------------------------------------------

#' Build the HuGeAMP BioIndex associations query URL for a given query string.
#'
#' @param query A variant rsID, varId (e.g. "1:100000:A:G"), or region
#'   (e.g. "1:100000-200000").
#' @param limit Maximum records to request (default 1000).
#' @return A fully-formed URL string.
build_hugeamp_query_url <- function(query, limit = 1000L) {
  sprintf(
    "%s/query/associations?q=%s&limit=%d",
    HUGEAMP_BASE,
    utils::URLencode(query, reserved = TRUE),
    as.integer(limit)
  )
}

#' Build the GWAS Catalog v2 associations query URL for a given variant.
#'
#' @param variant_id rsID or variant ID.
#' @param page_size  Results per page (default 500).
#' @return A fully-formed URL string.
build_gwas_query_url <- function(variant_id, page_size = 500L) {
  sprintf(
    "%s/associations?filter=variantId:%s&size=%d",
    GWAS_BASE,
    utils::URLencode(variant_id, reserved = FALSE),
    as.integer(page_size)
  )
}

# ---- Query diagnostics --------------------------------------------------------

#' Test a single API query URL and return a one-row summary data.frame.
#'
#' @param label      Short description of the query strategy being tested.
#' @param url        Fully-formed query URL.
#' @param parse_fn   A function that takes the response body (string) and the
#'   variant_id and returns a data.frame or NULL.
#' @param variant_id Used as a label / fallback in \code{parse_fn}.
#' @return A one-row data.frame with columns: strategy, url, http_status,
#'   n_results, error.
.test_query_url <- function(label, url, parse_fn, variant_id) {
  result_row <- data.frame(
    strategy    = label,
    url         = url,
    http_status = NA_integer_,
    n_results   = NA_integer_,
    error       = NA_character_,
    stringsAsFactors = FALSE
  )

  tryCatch({
    resp <- httr::GET(url,
                      httr::timeout(30),
                      httr::add_headers(Accept = "application/json"))

    result_row$http_status <- httr::status_code(resp)

    if (httr::http_error(resp)) {
      result_row$error <- paste0("HTTP ", httr::status_code(resp))
      return(result_row)
    }

    body <- httr::content(resp, "text", encoding = "UTF-8")
    df   <- tryCatch(parse_fn(body, variant_id), error = function(e) NULL)
    result_row$n_results <- if (is.null(df)) 0L else nrow(df)

  }, error = function(e) {
    result_row$error <<- conditionMessage(e)
  })

  result_row
}

#' Run a full set of query diagnostics for a variant across both APIs.
#'
#' Tests multiple query strategies in sequence and returns a data.frame
#' summarising the URL tried, the HTTP status code returned, the number of
#' association records found, and any error message.
#'
#' @param variant_id rsID or variant ID to test (e.g. "rs979012").
#' @return A data.frame with columns: strategy, url, http_status, n_results,
#'   error.  One row per strategy tested.
query_diagnostics <- function(variant_id) {
  rows <- list()

  # ---- GWAS Catalog -----------------------------------------------------------
  gwas_url <- build_gwas_query_url(variant_id)
  rows[["gwas_rsid"]] <- .test_query_url(
    label      = "GWAS Catalog: rsID query",
    url        = gwas_url,
    parse_fn   = function(body, vid) {
      data  <- jsonlite::fromJSON(body, simplifyDataFrame = FALSE)
      assoc <- data[["_embedded"]][["associations"]]
      if (is.null(assoc) || length(assoc) == 0) return(NULL)
      do.call(rbind, lapply(assoc, function(a) {
        data.frame(
          variant_id = safe_chr(a[["variant"]][["variantId"]], default = vid),
          study_id   = safe_chr(a[["study"]][["accessionId"]]),
          beta       = safe_num(a[["beta"]]),
          pvalue     = safe_num(a[["pValue"]]),
          stringsAsFactors = FALSE
        )
      }))
    },
    variant_id = variant_id
  )

  # ---- HuGeAMP: rsID query ----------------------------------------------------
  ha_rsid_url <- build_hugeamp_query_url(variant_id)
  rows[["hugeamp_rsid"]] <- .test_query_url(
    label      = "HuGeAMP: rsID query",
    url        = ha_rsid_url,
    parse_fn   = function(body, vid) .parse_hugeamp_body(body, vid),
    variant_id = variant_id
  )

  # ---- HuGeAMP: positional fallback via Ensembl ------------------------------
  is_rsid <- grepl(RSID_PATTERN, variant_id, ignore.case = TRUE)

  if (is_rsid &&
      !is.na(rows[["hugeamp_rsid"]]$http_status) &&
      rows[["hugeamp_rsid"]]$http_status == 400L) {

    coords <- resolve_rsid_to_grch38(variant_id)

    if (!is.null(coords)) {
      pos_query <- sprintf("%s:%d-%d", coords$chrom, coords$start, coords$end)
      ha_pos_url <- build_hugeamp_query_url(pos_query)
      rows[["hugeamp_position"]] <- .test_query_url(
        label      = paste0("HuGeAMP: positional query (", pos_query, ")"),
        url        = ha_pos_url,
        parse_fn   = function(body, vid) .parse_hugeamp_body(body, vid),
        variant_id = variant_id
      )
    } else {
      rows[["hugeamp_position"]] <- data.frame(
        strategy    = "HuGeAMP: positional query",
        url         = NA_character_,
        http_status = NA_integer_,
        n_results   = NA_integer_,
        error       = "Could not resolve rsID to GRCh38 position via Ensembl",
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
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
