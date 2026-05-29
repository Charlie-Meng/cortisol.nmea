#' Summarize benchmark winners
#'
#' @param summary_df Method-level benchmark summary data frame.
#' @param metrics Numeric metric columns where smaller is better.
#' @param group_vars Columns defining comparison groups.
#' @param method_col Method column name.
#' @param tol Tie tolerance.
#'
#' @return A data frame with one winner column per metric.
#' @export
summarize_winners <- function(summary_df,
                              metrics = c("Curve_RMSE", "AUC_MAE"),
                              group_vars = c("Scenario"),
                              method_col = "Method",
                              tol = 1e-12) {
  summary_df <- as.data.frame(summary_df)
  required <- c(group_vars, method_col, metrics)
  if (!all(required %in% names(summary_df))) {
    stop("`summary_df` is missing required columns.", call. = FALSE)
  }

  split_key <- interaction(summary_df[group_vars], drop = TRUE, lex.order = TRUE)
  groups <- split(summary_df, split_key)

  rows <- lapply(groups, function(dat) {
    out <- dat[1, group_vars, drop = FALSE]
    for (metric in metrics) {
      x <- dat[[metric]]
      ok <- is.finite(x)
      winner <- if (!any(ok)) {
        NA_character_
      } else {
        min_val <- min(x[ok], na.rm = TRUE)
        paste(dat[[method_col]][ok][abs(x[ok] - min_val) <= tol], collapse = ";")
      }
      out[[paste0("Winner_", metric)]] <- winner
    }
    out
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Build long and wide benchmark tables
#'
#' @param summary_df Method-level benchmark summary data frame.
#' @param metrics Numeric metric columns to keep.
#' @param id_cols Scenario identifier columns to keep.
#' @param method_col Method column name.
#'
#' @return A list with `long`, `wide`, and `winners` data frames.
#' @export
make_benchmark_tables <- function(summary_df,
                                  metrics = c("Curve_RMSE", "AUC_MAE"),
                                  id_cols = c("Scenario", "ScenarioGroup"),
                                  method_col = "Method") {
  summary_df <- as.data.frame(summary_df)
  id_cols <- intersect(id_cols, names(summary_df))
  required <- c(id_cols, method_col, metrics)
  if (!all(required %in% names(summary_df))) {
    stop("`summary_df` is missing required columns.", call. = FALSE)
  }

  long <- summary_df[, c(id_cols, method_col, metrics), drop = FALSE]
  metric_long <- do.call(rbind, lapply(metrics, function(metric) {
    data.frame(
      long[, c(id_cols, method_col), drop = FALSE],
      Metric = metric,
      Value = long[[metric]],
      row.names = NULL
    )
  }))

  wide <- stats::reshape(
    metric_long,
    idvar = c(id_cols, "Metric"),
    timevar = method_col,
    direction = "wide"
  )
  rownames(wide) <- NULL

  list(
    long = metric_long,
    wide = wide,
    winners = summarize_winners(summary_df, metrics = metrics, group_vars = id_cols, method_col = method_col)
  )
}
