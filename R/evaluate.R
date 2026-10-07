#' Evaluate a simulation study
#'
#' Summaries of [sim_run()] results. All errors compare fitted and true curves
#' on the design's evaluation grid (default 0 to 18 h in steps of 0.1 h).
#'
#' * `eval_subjects()`: one row per case, method and subject with the outcome,
#'   curve RMSE (nmol/L) and absolute error of the curve integral (curve-AUC
#'   error, nmol*h/L, trapezoidal rule). Subjects without a finite curve are
#'   unavailable and are never counted as recovered.
#' * `eval_summary()`: within each cohort, average subject errors over the
#'   subjects available for **all** compared methods, then summarize cohorts
#'   (mean, Monte Carlo SE, 95% t interval) by trimester and scenario, plus an
#'   equal-weight pooled row across trimesters (Welch degrees of freedom). The
#'   pooled row assumes independent cohorts across trimesters (see the `seeds`
#'   argument of [sim_design()]) and is `NA` when a planned trimester has no
#'   valid cohort, rather than silently re-weighting the others. A cohort in which any compared method failed is excluded from that
#'   comparison and counted in `Valid` vs `Planned`.
#' * `eval_paired()`: paired cohort-level differences `method - comparator` on
#'   the subjects common to the two methods; `Wins` counts cohorts where
#'   `method` has the lower error.
#' * `eval_detection()`: point-level detection of contaminated measurements by
#'   a workflow's outlier step (FPR, FNR, precision, and ROC-AUC of the
#'   absolute initial residual), per cohort and summarized.
#' * `eval_retention()`: percentage of subjects retained or excluded at each
#'   workflow step, per cohort and averaged.
#'
#' @param runs An `nmea_sim_runs` object.
#' @param methods Methods to compare (default: all).
#' @param method,comparator Method names.
#' @param metric `"rmse"` or `"auc_error"`.
#'
#' @return Data frames; `eval_summary()`, `eval_paired()`, `eval_detection()`
#'   and `eval_retention()` return a list with `cohorts` and `summary`.
#' @name eval_simulation
NULL

#' @rdname eval_simulation
#' @export
eval_subjects <- function(runs) {
  .check_runs(runs)
  times <- runs$design$times
  do.call(rbind, lapply(runs$cases, function(cs) {
    do.call(rbind, lapply(names(cs$methods), function(m) {
      r <- cs$methods[[m]]
      err <- r$curves - cs$truth_curves
      ok <- apply(r$curves, 1, function(x) all(is.finite(x)))
      data.frame(case = cs$id, trimester = cs$trimester, seed = cs$seed, scenario = cs$scenario,
                 method = m, subject = rownames(cs$truth_curves),
                 outcome = unname(r$outcome), available = unname(ok),
                 rmse = ifelse(ok, sqrt(rowMeans(err^2)), NA_real_),
                 auc_error = ifelse(ok, apply(err, 1, function(x) abs(.trapezoid(times, x))), NA_real_),
                 stringsAsFactors = FALSE)
    }))
  }))
}

#' @rdname eval_simulation
#' @export
eval_summary <- function(runs, methods = NULL) {
  .check_runs(runs)
  methods <- .check_methods(runs, methods)
  subj <- eval_subjects(runs)
  subj <- subj[subj$method %in% methods, ]
  cohorts <- do.call(rbind, lapply(runs$cases, function(cs) {
    failed <- any(vapply(cs$methods[methods], function(r) r$status != "ok", logical(1)))
    s <- subj[subj$case == cs$id, ]
    common <- Reduce(intersect, lapply(methods, function(m) s$subject[s$method == m & s$available]))
    do.call(rbind, lapply(methods, function(m) {
      q <- s[s$method == m & s$subject %in% common, ]
      data.frame(case = cs$id, trimester = cs$trimester, seed = cs$seed, scenario = cs$scenario,
                 method = m, valid = !failed, common_n = length(common),
                 own_n = sum(s$method == m & s$available),
                 rmse = if (!failed && nrow(q)) mean(q$rmse) else NA_real_,
                 auc_error = if (!failed && nrow(q)) mean(q$auc_error) else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }))
  summary <- .summarize_cohorts(cohorts, c("rmse", "auc_error", "common_n"), "method", valid = "valid")
  list(cohorts = cohorts, summary = summary)
}

#' @rdname eval_simulation
#' @export
eval_paired <- function(runs, method, comparator, metric = c("rmse", "auc_error")) {
  .check_runs(runs)
  metric <- match.arg(metric)
  .check_methods(runs, c(method, comparator))
  subj <- eval_subjects(runs)
  cohorts <- do.call(rbind, lapply(runs$cases, function(cs) {
    a <- subj[subj$case == cs$id & subj$method == method, ]
    b <- subj[subj$case == cs$id & subj$method == comparator, ]
    ok <- a$available & b$available[match(a$subject, b$subject)]
    valid <- cs$methods[[method]]$status == "ok" && cs$methods[[comparator]]$status == "ok" && any(ok)
    data.frame(case = cs$id, trimester = cs$trimester, seed = cs$seed, scenario = cs$scenario,
               valid = valid, common_n = sum(ok),
               difference = if (valid) mean(a[[metric]][ok] - b[[metric]][match(a$subject[ok], b$subject)]) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  cohorts$method <- paste(method, "-", comparator)
  summary <- .summarize_cohorts(cohorts, "difference", "method", valid = "valid", wins = TRUE)
  list(cohorts = cohorts, summary = summary)
}

#' @rdname eval_simulation
#' @export
eval_detection <- function(runs, method = "Full") {
  .check_runs(runs)
  .check_outlier_method(runs, method)
  cohorts <- do.call(rbind, lapply(runs$cases, function(cs) {
    r <- cs$methods[[method]]
    if (all(is.na(r$points$initial_residual))) {
      # The initial fit failed: no detection decision exists for this cohort.
      return(data.frame(case = cs$id, trimester = cs$trimester, seed = cs$seed, scenario = cs$scenario,
                        method = method, valid = FALSE, TP = NA, FP = NA, FN = NA, TN = NA,
                        cutoff = NA, FPR = NA, FNR = NA, precision = NA, roc_auc = NA))
    }
    truth <- cs$points$contaminated[match(r$points$point_id, cs$points$point_id)]
    f <- r$points$flagged
    tp <- sum(f & truth); fp <- sum(f & !truth); fn <- sum(!f & truth); tn <- sum(!f & !truth)
    data.frame(case = cs$id, trimester = cs$trimester, seed = cs$seed, scenario = cs$scenario,
               method = method, valid = TRUE, TP = tp, FP = fp, FN = fn, TN = tn, cutoff = r$cutoff,
               FPR = if (fp + tn) fp / (fp + tn) else NA_real_,
               FNR = if (tp + fn) fn / (tp + fn) else NA_real_,
               precision = if (tp + fp) tp / (tp + fp) else NA_real_,
               roc_auc = .roc_auc(abs(r$points$initial_residual), truth),
               stringsAsFactors = FALSE)
  }))
  summary <- .summarize_cohorts(cohorts, c("FPR", "FNR", "precision", "roc_auc", "cutoff"), "method",
                                valid = "valid")
  list(cohorts = cohorts, summary = summary)
}

#' @rdname eval_simulation
#' @export
eval_retention <- function(runs, method = "Full") {
  .check_runs(runs)
  .check_methods(runs, method)
  levels <- c("retained", "excluded_min_obs", "excluded_fvu", "excluded_c1", "not_fitted")
  cohorts <- do.call(rbind, lapply(runs$cases, function(cs) {
    r <- cs$methods[[method]]
    outcome <- unname(r$outcome)
    unavailable <- outcome == "retained" & !apply(r$curves, 1, function(x) all(is.finite(x)))
    outcome[unavailable] <- "not_fitted"
    pct <- 100 * as.numeric(table(factor(outcome, levels))) / length(outcome)
    cbind(data.frame(case = cs$id, trimester = cs$trimester, seed = cs$seed, scenario = cs$scenario,
                     method = method, valid = TRUE, stringsAsFactors = FALSE),
          stats::setNames(as.data.frame(t(pct)), levels))
  }))
  summary <- .summarize_cohorts(cohorts, levels, "method", valid = "valid")
  list(cohorts = cohorts, summary = summary)
}

# ---- internal helpers -------------------------------------------------------

.check_runs <- function(runs) {
  if (!inherits(runs, "nmea_sim_runs")) stop("`runs` must come from sim_run().", call. = FALSE)
}

.check_outlier_method <- function(runs, method) {
  .check_methods(runs, method)
  m <- runs$design$methods[[method]]
  if (m$type != "nmea" || is.null(m$steps$outlier)) {
    stop("Method `", method, "` has no outlier step.", call. = FALSE)
  }
  invisible(method)
}

.check_methods <- function(runs, methods) {
  available <- names(runs$design$methods)
  if (is.null(methods)) return(available)
  bad <- setdiff(methods, available)
  if (length(bad)) stop("Unknown method(s): ", paste(bad, collapse = ", "), call. = FALSE)
  methods
}

.trapezoid <- function(t, y) sum(diff(t) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

.roc_auc <- function(score, truth) {
  ok <- is.finite(score) & !is.na(truth)
  s <- score[ok]; t <- truth[ok]
  n1 <- sum(t); n0 <- sum(!t)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(s, ties.method = "average")[t]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Mean, Monte Carlo SE and 95% t interval of cohort values.
.mc <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  se <- if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
  q <- if (n > 1) stats::qt(0.975, n - 1) else NA_real_
  c(cohorts = n, mean = m, mcse = se, lower = m - q * se, upper = m + q * se)
}

# Equal-weight pooled mean across trimesters with a Welch-Satterthwaite interval.
# If a planned trimester has no valid cohort, the equal-weight target is not
# estimable: return NA instead of re-weighting the remaining trimesters.
.pooled <- function(x, trimester) {
  ok <- is.finite(x)
  groups <- split(x[ok], factor(trimester[ok], levels = sort(unique(trimester))))
  if (!length(groups) || any(lengths(groups) == 0)) {
    return(c(cohorts = sum(ok), mean = NA, mcse = NA, lower = NA, upper = NA))
  }
  n <- lengths(groups)
  w <- 1 / length(groups)
  m <- sum(w * vapply(groups, mean, numeric(1)))
  v <- w^2 * vapply(groups, function(g) if (length(g) > 1) stats::var(g) else NA_real_, numeric(1)) / n
  se <- sqrt(sum(v))
  df <- if (is.finite(se) && se > 0) sum(v)^2 / sum(v^2 / (n - 1)) else Inf
  q <- stats::qt(0.975, df)
  c(cohorts = sum(n), mean = m, mcse = se, lower = m - q * se, upper = m + q * se)
}

.summarize_cohorts <- function(cohorts, values, by, valid, wins = FALSE) {
  keys <- unique(cohorts[, c("scenario", by), drop = FALSE])
  rows <- list()
  for (i in seq_len(nrow(keys))) {
    sel <- cohorts$scenario == keys$scenario[i] & cohorts[[by]] == keys[[by]][i]
    d <- cohorts[sel, ]
    for (v in values) {
      x <- ifelse(d[[valid]], d[[v]], NA_real_)
      for (tri in sort(unique(d$trimester))) {
        s <- .mc(x[d$trimester == tri])
        rows[[length(rows) + 1]] <- data.frame(trimester = tri, scenario = keys$scenario[i],
          method = keys[[by]][i], measure = v, planned = sum(d$trimester == tri), t(s),
          wins = if (wins) sum(x[d$trimester == tri] < 0, na.rm = TRUE) else NA_integer_)
      }
      if (length(unique(d$trimester)) > 1) {
        s <- .pooled(x, d$trimester)
        rows[[length(rows) + 1]] <- data.frame(trimester = "All", scenario = keys$scenario[i],
          method = keys[[by]][i], measure = v, planned = nrow(d), t(s),
          wins = if (wins) sum(x < 0, na.rm = TRUE) else NA_integer_)
      }
    }
  }
  out <- do.call(rbind, rows)
  names(out)[names(out) == "cohorts"] <- "valid"
  rownames(out) <- NULL
  out
}
