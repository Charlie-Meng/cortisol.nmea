#' Fit a symmetric baseline NMEA model
#'
#' Convenience wrapper around [fit_nmea()] with `alpha = 1`.
#'
#' @inheritParams fit_nmea
#'
#' @return An `nmea_fit` object.
#' @export
fit_baseline_nmea <- function(data,
                              id = "Subject",
                              time = "Time",
                              response = "Conc_Obs",
                              psi0 = c(mu = 1, s = 1, c1 = 1, c0 = 0),
                              seed = 632545,
                              quiet = TRUE,
                              ...) {
  fit_nmea(
    data = data,
    id = id,
    time = time,
    response = response,
    alpha = 1,
    psi0 = psi0,
    seed = seed,
    quiet = quiet,
    ...
  )
}

#' Fit an alpha-grid NMEA model
#'
#' Convenience wrapper around [fit_alpha_grid()] that returns the selected final
#' fit with the alpha-grid object attached as `alpha_grid_result`.
#'
#' @inheritParams fit_alpha_grid
#'
#' @return An `nmea_fit` object.
#' @export
fit_alpha_only_nmea <- function(data,
                                alpha_grid = seq(0.8, 1.5, by = 0.1),
                                id = "Subject",
                                time = "Time",
                                response = "Conc_Obs",
                                psi0 = c(mu = 1, s = 1, c1 = 1, c0 = 0),
                                seed = 632545,
                                quiet = TRUE,
                                ...) {
  grid_fit <- fit_alpha_grid(
    data = data,
    alpha_grid = alpha_grid,
    id = id,
    time = time,
    response = response,
    psi0 = psi0,
    seed = seed,
    quiet = quiet,
    ...
  )
  out <- grid_fit$final_fit
  out$alpha_grid_result <- grid_fit
  out
}

#' Fit a smooth GAM baseline
#'
#' Fits a lightweight smooth competitor using `mgcv::gam()`. This is a public
#' package baseline for synthetic benchmarks, not a full replacement for the
#' project-specific Brisa-style GAMM analysis scripts.
#'
#' @param data Data frame with subject, time, and response columns.
#' @param id Subject identifier column.
#' @param time Time column.
#' @param response Response column.
#' @param k Basis size for the population smooth.
#' @param subject_effect Include a random subject intercept with `s(subject,
#'   bs = "re")`.
#' @param method Smoothing estimation method passed to `mgcv::gam()`.
#'
#' @return A list with the fitted model and column metadata.
#' @export
fit_gam_baseline <- function(data,
                             id = "Subject",
                             time = "Time",
                             response = "Conc_Obs",
                             k = 5,
                             subject_effect = TRUE,
                             method = "REML") {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package `mgcv` is required for `fit_gam_baseline()`.", call. = FALSE)
  }
  data <- as.data.frame(data)
  if (!all(c(id, time, response) %in% names(data))) {
    stop("`data` must contain id, time, and response columns.", call. = FALSE)
  }

  dat <- data.frame(
    .subject = factor(data[[id]]),
    .time = as.numeric(data[[time]]),
    .response = as.numeric(data[[response]])
  )

  formula <- if (isTRUE(subject_effect)) {
    stats::as.formula(paste0(".response ~ s(.time, k = ", k, ") + s(.subject, bs = 're')"))
  } else {
    stats::as.formula(paste0(".response ~ s(.time, k = ", k, ")"))
  }

  fit <- mgcv::gam(formula, data = dat, method = method)
  structure(
    list(
      fit = fit,
      training_subjects = levels(dat$.subject),
      columns = list(id = id, time = time, response = response),
      subject_effect = subject_effect
    ),
    class = "nmea_gam_baseline"
  )
}

#' Predict from a smooth GAM baseline
#'
#' @param object Object returned by [fit_gam_baseline()].
#' @param subject_ids Character vector of subject IDs to predict for.
#' @param time_grid Numeric vector of prediction times.
#'
#' @return Numeric prediction matrix.
#' @export
predict_gam_baseline <- function(object, subject_ids, time_grid) {
  if (!inherits(object, "nmea_gam_baseline")) {
    stop("`object` must be returned by `fit_gam_baseline()`.", call. = FALSE)
  }
  subject_ids <- as.character(subject_ids)
  newdata <- expand.grid(
    .time = time_grid,
    .subject = subject_ids,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  newdata$.subject <- factor(newdata$.subject, levels = object$training_subjects)
  pred <- stats::predict(object$fit, newdata = newdata)
  out <- matrix(pred, nrow = length(time_grid), ncol = length(subject_ids))
  out <- t(out)
  rownames(out) <- subject_ids
  colnames(out) <- paste0("t=", time_grid)
  out
}
