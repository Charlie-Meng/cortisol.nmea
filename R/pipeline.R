#' Steps of the NMEA workflow
#'
#' Configures which steps [nmea_pipeline()] runs. Setting a step to `NULL`
#' (or `FALSE` for `c1_positive`) skips it.
#'
#' The full workflow is:
#' 1. initial fit with `alpha = initial_alpha`;
#' 2. remove measurements flagged by the `outlier` rule;
#' 3. exclude subjects left with fewer than `min_obs` measurements;
#' 4. refit at `initial_alpha` and exclude subjects whose fraction of variance
#'    unexplained is not below `fvu_max`, and (if `c1_positive`) subjects
#'    whose amplitude `c1` is not positive, i.e. curves without a morning peak;
#' 5. fit over `alpha_grid` and keep the alpha with the smallest residual sum
#'    of squares (a single value fixes alpha).
#'
#' The refit in step 4 runs only when at least one screen is enabled.
#'
#' @param initial_alpha Alpha for the initial and screening fits.
#' @param outlier An outlier rule (see [outlier_fixed()]), a number (treated
#'   as a fixed cutoff), or `NULL` to keep all measurements.
#' @param min_obs Minimum number of measurements per subject after outlier
#'   removal, or `NULL`.
#' @param fvu_max Subjects need a fraction of variance unexplained below this
#'   value, or `NULL`.
#' @param c1_positive Logical; exclude subjects with non-positive amplitude.
#' @param alpha_grid Candidate alpha values for the final fit.
#'
#' @return An object of class `nmea_steps`.
#' @export
#' @examples
#' nmea_steps()
#' nmea_steps(outlier = outlier_sd(k = 3.5, sd = "mad"), fvu_max = NULL)
#' nmea_steps_direct()
nmea_steps <- function(initial_alpha = 1,
                       outlier = outlier_fixed(6),
                       min_obs = 3L,
                       fvu_max = 0.5,
                       c1_positive = TRUE,
                       alpha_grid = seq(0.8, 1.5, by = 0.1)) {
  .check_scalar_number(initial_alpha, "initial_alpha", lower = .Machine$double.eps)
  if (is.numeric(outlier)) outlier <- outlier_fixed(outlier)
  if (!is.null(outlier) && !inherits(outlier, "nmea_outlier_rule")) {
    stop("`outlier` must be an nmea_outlier_rule, a number, or NULL.", call. = FALSE)
  }
  if (!is.null(min_obs)) .check_scalar_number(min_obs, "min_obs", lower = 1)
  if (!is.null(fvu_max)) .check_scalar_number(fvu_max, "fvu_max", lower = 0)
  if (!is.logical(c1_positive) || length(c1_positive) != 1L || is.na(c1_positive)) {
    stop("`c1_positive` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(alpha_grid) || !length(alpha_grid) || any(!is.finite(alpha_grid)) ||
      any(alpha_grid <= 0)) {
    stop("`alpha_grid` must contain positive finite values.", call. = FALSE)
  }
  structure(list(initial_alpha = initial_alpha, outlier = outlier, min_obs = min_obs,
                 fvu_max = fvu_max, c1_positive = c1_positive, alpha_grid = alpha_grid),
            class = "nmea_steps")
}

#' @rdname nmea_steps
#' @export
nmea_steps_direct <- function(alpha_grid = seq(0.8, 1.5, by = 0.1)) {
  nmea_steps(outlier = NULL, min_obs = NULL, fvu_max = NULL, c1_positive = FALSE,
             alpha_grid = alpha_grid)
}

#' @export
print.nmea_steps <- function(x, ...) {
  cat("<nmea_steps>\n")
  cat("  initial alpha :", format(x$initial_alpha), "\n")
  cat("  outlier rule  :", if (is.null(x$outlier)) "none" else x$outlier$label, "\n")
  cat("  min obs       :", if (is.null(x$min_obs)) "none" else x$min_obs, "\n")
  cat("  FVU screen    :", if (is.null(x$fvu_max)) "none" else paste("<", x$fvu_max), "\n")
  cat("  c1 > 0 screen :", x$c1_positive, "\n")
  cat("  alpha grid    :", paste(format(x$alpha_grid), collapse = ", "), "\n")
  invisible(x)
}

#' Run the NMEA workflow
#'
#' Runs the configurable workflow described in [nmea_steps()] and records, for
#' every measurement and subject, what happened at each step. Measurement flags
#' and subject exclusions are kept separate: a subject can lose all of its
#' curve because too few measurements remained, even if most were not flagged.
#'
#' Fit each trimester (or other homogeneous group) separately.
#'
#' @param data A [cortisol_data()] object.
#' @param steps Workflow configuration from [nmea_steps()]; use
#'   [nmea_steps_direct()] for an uncleaned fit of all data.
#' @param control Fitting settings from [nmea_control()].
#' @param keep_fits Logical; keep every alpha-grid fit.
#'
#' @return An object of class `nmea_result` with elements:
#'   * `status`: `"ok"` or `"failure"`, and `error`;
#'   * `final`: the selected final `nmea_fit`, and `alpha`, `profile`;
#'   * `points`: one row per input measurement with the initial residual,
#'     the flag and the cutoff;
#'   * `subjects`: one row per input subject with its outcome (`retained`,
#'     `excluded_min_obs`, `excluded_fvu`, `excluded_c1` or `not_fitted`);
#'   * `stages`: the initial and screening fits;
#'   * `steps`, `control`.
#' @export
nmea_pipeline <- function(data, steps = nmea_steps(), control = nmea_control(),
                          keep_fits = FALSE) {
  data <- .as_cortisol_data(data)
  if (!inherits(steps, "nmea_steps")) stop("`steps` must come from nmea_steps().", call. = FALSE)
  if ("trimester" %in% names(data) && length(unique(data$trimester)) > 1L) {
    warning("`data` contains several trimesters; NMEA is meant to be fitted to each separately.",
            call. = FALSE)
  }

  subjects <- sort(unique(data$subject))
  points <- data.frame(point_id = data$point_id, subject = data$subject, time = data$time,
                       y = data$y, initial_residual = NA_real_, flagged = FALSE,
                       stringsAsFactors = FALSE)
  subject_table <- data.frame(subject = subjects,
                              n_obs = as.integer(table(data$subject)[subjects]),
                              n_flagged = 0L, outcome = "not_fitted",
                              stringsAsFactors = FALSE)
  result <- list(status = "failure", error = NULL, final = NULL, alpha = NA_real_,
                 profile = NULL, points = points, subjects = subject_table,
                 cutoff = NA_real_, stages = list(), steps = steps, control = control)
  finish <- function(result, error = NULL) {
    result$error <- error
    if (!is.null(error)) result$status <- "failure"
    structure(result, class = "nmea_result")
  }

  current <- data
  # Step 1-2: initial fit and outlier removal.
  if (!is.null(steps$outlier)) {
    first <- nmea_fit(current, alpha = steps$initial_alpha, control = control)
    result$stages$initial <- first
    if (first$status != "ok") return(finish(result, paste("Initial fit failed:", first$error)))
    rule <- .apply_outlier_rule(steps$outlier, first)
    idx <- match(result$points$point_id, first$data$point_id)
    result$points$initial_residual <- first$residual[idx]
    result$points$flagged <- rule$flag[idx]
    result$cutoff <- rule$cutoff
    flagged_by_subject <- tapply(result$points$flagged, result$points$subject, sum)
    result$subjects$n_flagged <- as.integer(flagged_by_subject[subjects])
    current <- .subset_cd(first$data, !rule$flag)
  }

  # Step 3: minimum number of measurements.
  if (!is.null(steps$min_obs)) {
    counts <- table(current$subject)
    eligible <- names(counts)[counts >= steps$min_obs]
    removed <- setdiff(subjects, eligible)
    result$subjects$outcome[result$subjects$subject %in% removed] <- "excluded_min_obs"
    current <- .subset_cd(current, current$subject %in% eligible)
  }
  if (length(unique(current$subject)) < 3L) {
    return(finish(result, "Fewer than 3 subjects remain after measurement removal."))
  }

  # Step 4: screening refit.
  if (!is.null(steps$fvu_max) || isTRUE(steps$c1_positive)) {
    screen <- nmea_fit(current, alpha = steps$initial_alpha, control = control)
    result$stages$screening <- screen
    if (screen$status != "ok") return(finish(result, paste("Screening fit failed:", screen$error)))
    keep <- rownames(screen$psi)
    if (!is.null(steps$fvu_max)) {
      fvu_ok <- names(screen$fvu)[is.finite(screen$fvu) & screen$fvu < steps$fvu_max]
      result$subjects$outcome[result$subjects$subject %in% setdiff(keep, fvu_ok)] <- "excluded_fvu"
      keep <- intersect(keep, fvu_ok)
    }
    if (isTRUE(steps$c1_positive)) {
      c1_ok <- rownames(screen$psi)[screen$psi[, "c1"] > 0]
      result$subjects$outcome[result$subjects$subject %in% setdiff(keep, c1_ok)] <- "excluded_c1"
      keep <- intersect(keep, c1_ok)
    }
    current <- .subset_cd(current, current$subject %in% keep)
    if (length(keep) < 3L) {
      return(finish(result, "Fewer than 3 subjects remain after screening."))
    }
  }

  # Step 5: final fit over the alpha grid.
  final <- nmea_fit_alpha(current, alpha_grid = steps$alpha_grid, control = control,
                          keep_fits = keep_fits)
  result$profile <- final$profile
  if (keep_fits) result$stages$alpha_fits <- final$fits
  if (final$status != "ok") return(finish(result, paste("Final fit:", final$error)))
  result$status <- "ok"
  result$final <- final$best
  result$alpha <- final$alpha
  result$subjects$outcome[result$subjects$subject %in% rownames(final$best$psi)] <- "retained"
  finish(result)
}

#' @export
print.nmea_result <- function(x, ...) {
  cat("<nmea_result>", if (x$status == "ok") "completed" else paste("FAILED:", x$error), "\n")
  cat(sprintf("  %d subjects, %d measurements; %d measurements flagged",
              nrow(x$subjects), nrow(x$points), sum(x$points$flagged)))
  if (is.finite(x$cutoff)) cat(sprintf(" (cutoff %.3f nmol/L)", x$cutoff))
  cat("\n")
  tab <- table(factor(x$subjects$outcome, levels = c("retained", "excluded_min_obs",
                                                     "excluded_fvu", "excluded_c1",
                                                     "not_fitted")))
  cat("  subject outcomes:", paste(names(tab), tab, sep = " = ", collapse = ", "), "\n")
  if (x$status == "ok") cat("  selected alpha:", format(x$alpha), "\n")
  invisible(x)
}

#' @export
coef.nmea_result <- function(object, ...) {
  if (object$status != "ok") stop("The workflow failed; no parameters are available.", call. = FALSE)
  object$final$psi
}

#' Predict NMEA curves
#'
#' Evaluates fitted curves on a time grid with the numerically stable
#' [nmea_curve()].
#'
#' @param object An `nmea_result` or `nmea_fit`.
#' @param times Times since waking (hours).
#' @param subjects Subject identifiers to return, always matched by identifier
#'   (numeric identifiers are matched as text, never used as row positions).
#'   For an `nmea_result` the default is every input subject, and subjects
#'   without a final curve get rows of `NA`; unknown identifiers are an error.
#'   For an `nmea_fit` the default is every fitted subject.
#' @param ... Unused.
#'
#' @return A matrix with one row per subject and one column per time.
#' @export
predict.nmea_result <- function(object, times = seq(0, 18, by = 0.1), subjects = NULL, ...) {
  subjects <- if (is.null(subjects)) object$subjects$subject else as.character(subjects)
  .check_known_subjects(subjects, object$subjects$subject)
  out <- matrix(NA_real_, length(subjects), length(times), dimnames = list(subjects, NULL))
  if (object$status != "ok") return(out)
  psi <- object$final$psi
  have <- intersect(subjects, rownames(psi))
  if (length(have)) out[have, ] <- nmea_curve(times, psi[have, , drop = FALSE])
  out
}

#' @rdname predict.nmea_result
#' @export
predict.nmea_fit <- function(object, times = seq(0, 18, by = 0.1), subjects = NULL, ...) {
  psi <- stats::coef(object)
  subjects <- if (is.null(subjects)) rownames(psi) else as.character(subjects)
  .check_known_subjects(subjects, rownames(psi))
  nmea_curve(times, psi[match(subjects, rownames(psi)), , drop = FALSE])
}

.check_known_subjects <- function(subjects, known) {
  unknown <- setdiff(subjects, known)
  if (length(unknown)) {
    stop("Unknown subject identifier(s): ", paste(utils::head(unknown, 5), collapse = ", "),
         if (length(unknown) > 5) ", ...", call. = FALSE)
  }
  invisible(subjects)
}
