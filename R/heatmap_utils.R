# R/heatmap_utils.R
# Utilities for constructing ComplexHeatmap objects from beta-value matrices.

#' Build a ComplexHeatmap::Heatmap from a variants × studies beta matrix.
#'
#' Rows    = variants (rsIDs)
#' Columns = study identifiers
#' Values  = association beta values (NA rendered as grey)
#'
#' @param beta_matrix Numeric matrix with rownames (variants) and colnames
#'   (studies).
#' @return A \code{\link[ComplexHeatmap]{Heatmap}} object, or NULL if the
#'   matrix is empty.
build_heatmap <- function(beta_matrix) {
  if (is.null(beta_matrix) ||
      nrow(beta_matrix) == 0 ||
      ncol(beta_matrix) == 0) {
    return(NULL)
  }

  # Determine symmetric colour scale centred on 0
  finite_vals   <- beta_matrix[is.finite(beta_matrix)]
  max_abs_beta  <- if (length(finite_vals) > 0) max(abs(finite_vals)) else 1
  if (max_abs_beta == 0) max_abs_beta <- 1

  col_fun <- circlize::colorRamp2(
    c(-max_abs_beta, 0, max_abs_beta),
    c("#3182BD", "white", "#E6550D")   # blue → white → orange
  )

  # Decide whether to cluster (need ≥ 2 rows/cols with non-NA data)
  do_row_cluster <- nrow(beta_matrix) > 2 &&
                    sum(rowSums(!is.na(beta_matrix)) > 0) > 2
  do_col_cluster <- ncol(beta_matrix) > 2 &&
                    sum(colSums(!is.na(beta_matrix)) > 0) > 2

  ComplexHeatmap::Heatmap(
    beta_matrix,
    name                  = "Beta",
    col                   = col_fun,
    na_col                = "grey85",
    rect_gp               = grid::gpar(col = "white", lwd = 0.8),
    row_title             = "Variants",
    column_title          = "Studies",
    row_names_gp          = grid::gpar(fontsize = 10),
    column_names_gp       = grid::gpar(fontsize = 9),
    column_names_rot      = 45,
    cluster_rows          = do_row_cluster,
    cluster_columns       = do_col_cluster,
    show_row_dend         = do_row_cluster,
    show_column_dend      = do_col_cluster,
    heatmap_legend_param  = list(
      title              = "Beta",
      title_gp           = grid::gpar(fontsize = 10, fontface = "bold"),
      direction          = "horizontal",
      legend_width       = grid::unit(4, "cm")
    )
  )
}
