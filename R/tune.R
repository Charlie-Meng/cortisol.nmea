#' Threshold trade-offs without refitting
#'
#' Re-applies `cutoff = k x sd` to the saved initial residuals of a workflow
#' method and reports, per cohort and on average, the false positive rate
#' (clean measurements flagged) and false negative rate (contaminated
#' measurements missed) for every `k`. Because the initial fit does not depend
#' on the cutoff, no model is refitted; downstream curve recovery and subject
#' retention are **not** evaluated. To compare thresholds including refitting,
#' run a design with [methods_threshold_grid()].
#'
#' The function reports trade-offs only; it does not choose a threshold.
#'
#' @param runs An `nmea_sim_runs` object.
#' @param k Numeric vector of multipliers.
#' @param sd Noise scale: a number, `"reference"` (the trimester reference
#'   noise SD of the design), `"model"` (initial-fit residual SD) or `"mad"`
#'   (MAD of initial residuals).
#' @param method Name of a workflow method with an outlier step.
#'
#' @return A list with `cohorts` (one row per case and `k`; cases whose initial
#'   fit failed are kept with `valid = FALSE` and the error) and `summary`
#'   (planned and valid cohorts, mean FPR and FNR by trimester, scenario and `k`).
#' @export
tune_threshold <- function(runs, k, sd = "reference", method = "Full") {
  .check_runs(runs)
  .check_outlier_method(runs, method)
  if (!is.numeric(k) || any(!is.finite(k)) || any(k <= 0)) stop("`k` must be positive.", call. = FALSE)
  cohorts <- do.call(rbind, lapply(runs$cases, function(cs) {
    r <- cs$methods[[method]]
    res <- r$points$initial_residual
    if (all(is.na(res))) {
      # Initial fit failed: keep the planned rows as invalid, with the reason.
      return(data.frame(case = cs$id, trimester = cs$trimester, seed = cs$seed, scenario = cs$scenario,
                        k = k, valid = FALSE, cutoff = NA_real_, FPR = NA_real_, FNR = NA_real_,
                        flagged = NA_integer_, error = if (is.null(r$error)) NA_character_ else r$error,
                        stringsAsFactors = FALSE))
    }
    truth <- cs$points$contaminated[match(r$points$point_id, cs$points$point_id)]
    scale <- if (is.numeric(sd)) sd else switch(sd,
      reference = runs$design$references[[cs$trimester]]$noise_sd,
      model = r$initial_sigma,
      mad = stats::mad(res),
      iterative = stop("sd = \"iterative\" needs refitting; compare it with methods_threshold_grid().",
                       call. = FALSE),
      stop("`sd` must be a number, \"reference\", \"model\" or \"mad\".", call. = FALSE))
    do.call(rbind, lapply(k, function(kk) {
      f <- abs(res) > kk * scale
      data.frame(case = cs$id, trimester = cs$trimester, seed = cs$seed, scenario = cs$scenario,
                 k = kk, valid = TRUE, cutoff = kk * scale,
                 FPR = if (any(!truth)) mean(f[!truth]) else NA_real_,
                 FNR = if (any(truth)) mean(!f[truth]) else NA_real_,
                 flagged = sum(f), error = NA_character_, stringsAsFactors = FALSE)
    }))
  }))
  keys <- unique(cohorts[, c("trimester", "scenario", "k")])
  mean_or_na <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  summary <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    d <- cohorts[cohorts$trimester == keys$trimester[i] & cohorts$scenario == keys$scenario[i] &
                   cohorts$k == keys$k[i], ]
    v <- d[d$valid, ]
    data.frame(keys[i, ], planned = nrow(d), valid = nrow(v), cutoff = mean_or_na(v$cutoff),
               FPR = mean_or_na(v$FPR), FNR = mean_or_na(v$FNR), flagged = mean_or_na(v$flagged))
  }))
  summary <- summary[order(summary$scenario, summary$trimester, summary$k), ]
  rownames(summary) <- NULL
  list(cohorts = cohorts, summary = summary)
}
