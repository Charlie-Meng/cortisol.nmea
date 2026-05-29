#' Curve-level benchmark metrics
#'
#' Computes subject-level curve RMSE and AUC error from estimated and true curve
#' matrices on a shared time grid.
#'
#' @param estimated Matrix of estimated curves.
#' @param truth Matrix of true curves.
#' @param time_grid Numeric time grid shared by both matrices.
#'
#' @return A data frame with subject-level benchmark metrics.
#' @export
curve_benchmark_metrics <- function(estimated, truth, time_grid) {
  common <- intersect(rownames(estimated), rownames(truth))
  if (length(common) == 0) {
    return(data.frame(
      Subject = character(0),
      Curve_RMSE = numeric(0),
      AUC_Error = numeric(0)
    ))
  }

  estimated <- estimated[common, , drop = FALSE]
  truth <- truth[common, , drop = FALSE]
  data.frame(
    Subject = common,
    Curve_RMSE = sqrt(rowMeans((estimated - truth)^2, na.rm = TRUE)),
    AUC_Error = vapply(seq_along(common), function(i) {
      .trapz(estimated[i, ], time_grid) - .trapz(truth[i, ], time_grid)
    }, numeric(1)),
    row.names = NULL
  )
}

#' Summarize subject-level benchmark metrics
#'
#' @param metric_df Output from `curve_benchmark_metrics()`.
#'
#' @return One-row data frame of aggregate metrics.
#' @export
summarize_curve_metrics <- function(metric_df) {
  data.frame(
    N_Subjects = nrow(metric_df),
    Curve_RMSE = mean(metric_df$Curve_RMSE, na.rm = TRUE),
    AUC_MAE = mean(abs(metric_df$AUC_Error), na.rm = TRUE)
  )
}
