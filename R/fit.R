#' Control settings for NMEA fitting
#'
#' @param iterations Length-two integer vector of SAEM iterations for the
#'   exploration and smoothing phases (saemix `nbiter.saemix`).
#' @param seed Seed passed to saemix. The global random number stream of the
#'   R session is restored after each fit.
#' @param psi0 Named starting values for `mu`, `s`, `c1` and `c0`.
#' @param fim Logical; compute the Fisher information matrix and standard
#'   errors.
#' @param quiet Logical; suppress saemix console output.
#'
#' @return A list of class `nmea_control`.
#' @export
#' @examples
#' nmea_control(iterations = c(150, 75))
nmea_control <- function(iterations = c(300L, 100L), seed = 632545L,
                         psi0 = c(mu = 1, s = 1, c1 = 1, c0 = 0),
                         fim = TRUE, quiet = TRUE) {
  if (!is.numeric(iterations) || length(iterations) != 2L || any(iterations < 1)) {
    stop("`iterations` must be two positive numbers.", call. = FALSE)
  }
  if (is.null(names(psi0)) || !setequal(names(psi0), c("mu", "s", "c1", "c0"))) {
    stop("`psi0` must be named mu, s, c1 and c0.", call. = FALSE)
  }
  structure(list(
    iterations = as.integer(iterations),
    seed = as.integer(seed),
    psi0 = psi0[c("mu", "s", "c1", "c0")],
    fim = isTRUE(fim),
    quiet = isTRUE(quiet)
  ), class = "nmea_control")
}

#' Fit the NMEA model with a fixed asymmetry parameter
#'
#' Fits the asymmetric scaled-logistic nonlinear mixed-effects model by SAEM
#' (package saemix) with diagonal random effects on `mu`, `s`, `c1` and `c0`,
#' a constant residual error model, and `alpha` held fixed. Individual
#' parameters are maximum a posteriori (MAP) estimates.
#'
#' Measurements are ordered by subject, time and measurement identifier before
#' fitting, so results do not depend on the input row order.
#'
#' @param data A [cortisol_data()] object.
#' @param alpha Fixed asymmetry parameter.
#' @param control Settings from [nmea_control()].
#'
#' @return An object of class `nmea_fit` with elements `status` (`"ok"` or
#'   `"failure"`), `alpha`, `psi` (individual parameters), `data` (ordered
#'   data used for fitting), `fitted`, `residual`, `rss`, `fvu` (per-subject
#'   fraction of variance unexplained), `sigma` (residual error SD), `fim_ok`,
#'   `saemix` (the saemix fit), `warnings`, `error` and `elapsed`.
#' @references Comets E, Lavenu A, Lavielle M (2017). Parameter estimation in
#'   nonlinear mixed effect models using saemix, an R implementation of the
#'   SAEM algorithm. *Journal of Statistical Software*, 80(3), 1-41.
#'   \doi{10.18637/jss.v080.i03}
#' @export
nmea_fit <- function(data, alpha = 1, control = nmea_control()) {
  data <- .as_cortisol_data(data)
  .check_scalar_number(alpha, "alpha", lower = .Machine$double.eps)
  if (!inherits(control, "nmea_control")) stop("`control` must come from nmea_control().", call. = FALSE)

  ids <- sort(unique(data$subject))
  ord <- order(match(data$subject, ids), data$time, data$point_id)
  dd <- .subset_cd(data, ord)
  started <- proc.time()[["elapsed"]]

  res <- .collect_conditions(.with_preserved_rng(.maybe_quiet(
    .nmea_fit_saemix(dd, ids, alpha, control), quiet = control$quiet
  )))
  elapsed <- proc.time()[["elapsed"]] - started

  if (!is.null(res$error)) {
    return(structure(list(status = "failure", error = res$error, alpha = alpha, data = dd,
                          warnings = res$warnings, elapsed = elapsed, control = control),
                     class = "nmea_fit"))
  }
  out <- res$value
  out$status <- "ok"
  out$error <- NULL
  out$warnings <- res$warnings
  out$elapsed <- elapsed
  out$control <- control
  structure(out, class = "nmea_fit")
}

.nmea_fit_saemix <- function(dd, ids, alpha, control) {
  group <- match(dd$subject, ids)
  sdat <- saemix::saemixData(
    name.data = data.frame(id = group, time = dd$time, y = dd$y),
    name.group = "id", name.predictors = "time", name.response = "y",
    units = list(x = "Hour", y = "nmol/L"), verbose = FALSE
  )
  model <- saemix::saemixModel(
    model = .nmea_saemix_model,
    psi0 = c(control$psi0, alpha = alpha),
    fixed.estim = c(1, 1, 1, 1, 0),
    covariance.model = diag(c(1, 1, 1, 1, 0)),
    verbose = FALSE
  )
  fit <- saemix::saemix(model, sdat, control = list(
    map = TRUE, fim = control$fim, ll.is = FALSE,
    nbiter.saemix = control$iterations, seed = control$seed,
    displayProgress = FALSE, print = FALSE, save = FALSE, save.graphs = FALSE
  ))

  psi <- as.matrix(saemix::psi(fit))
  index_map <- unique(fit@data@data[, c("index", "id")])
  index_map <- index_map[order(index_map$index), , drop = FALSE]
  if (nrow(psi) != nrow(index_map) || !setequal(index_map$id, seq_along(ids)) ||
      !identical(as.numeric(fit@data@data$time), as.numeric(dd$time)) ||
      !identical(as.numeric(fit@data@data$y), as.numeric(dd$y))) {
    stop("saemix reordered the data; subject mapping could not be verified.", call. = FALSE)
  }
  rownames(psi) <- ids[index_map$id]
  colnames(psi) <- .psi_names

  fitted <- as.numeric(fit@results@ipred)
  map_check <- .nmea_saemix_model(psi, match(dd$subject, rownames(psi)), matrix(dd$time, ncol = 1))
  if (length(fitted) != nrow(dd) || any(!is.finite(psi)) || any(!is.finite(fitted)) ||
      max(abs(fitted - map_check)) > 1e-6) {
    stop("Fitted values are not finite or do not match the MAP parameters.", call. = FALSE)
  }
  residual <- dd$y - fitted
  sse <- tapply(residual^2, dd$subject, sum)
  sst <- tapply(dd$y, dd$subject, function(v) sum((v - mean(v))^2))
  se <- fit@results@se.fixed

  list(
    alpha = alpha, psi = psi, data = dd, fitted = fitted, residual = residual,
    rss = sum(residual^2), fvu = sse / sst,
    sigma = as.numeric(fit@results@respar[1]),
    fim_ok = length(se) >= 4L && all(is.finite(se[1:4])) && all(se[1:4] > 0),
    saemix = fit
  )
}

#' @export
print.nmea_fit <- function(x, ...) {
  if (x$status != "ok") {
    cat(sprintf("<nmea_fit> FAILED at alpha = %s: %s\n", format(x$alpha), x$error))
    return(invisible(x))
  }
  cat(sprintf("<nmea_fit> alpha = %s, %d subjects, %d measurements\n",
              format(x$alpha), nrow(x$psi), nrow(x$data)))
  cat(sprintf("RSS = %.4f, residual SD = %.4f, FIM available: %s\n", x$rss, x$sigma, x$fim_ok))
  if (length(x$warnings)) cat(sprintf("%d warning(s) recorded in $warnings\n", length(x$warnings)))
  invisible(x)
}

#' @export
coef.nmea_fit <- function(object, ...) {
  if (object$status != "ok") stop("The fit failed; no parameters are available.", call. = FALSE)
  object$psi
}

#' Fit the NMEA model over a grid of asymmetry values
#'
#' Fits [nmea_fit()] for each `alpha` and selects the value with the smallest
#' residual sum of squares. Failed fits are recorded and never selected.
#'
#' @inheritParams nmea_fit
#' @param alpha_grid Numeric vector of candidate `alpha` values.
#' @param keep_fits Logical; keep every fit (memory intensive) or only the
#'   selected one.
#'
#' @return An object of class `nmea_alpha_fit` with `status`, `alpha`
#'   (selected), `best` (an `nmea_fit`), `profile` (one row per candidate with
#'   alpha, status, RSS, FIM availability, warnings, error and elapsed time,
#'   kept even when `keep_fits = FALSE`) and, optionally, `fits`.
#' @export
nmea_fit_alpha <- function(data, alpha_grid = seq(0.8, 1.5, by = 0.1),
                           control = nmea_control(), keep_fits = FALSE) {
  if (!is.numeric(alpha_grid) || !length(alpha_grid) || any(!is.finite(alpha_grid)) ||
      any(alpha_grid <= 0)) {
    stop("`alpha_grid` must contain positive finite values.", call. = FALSE)
  }
  data <- .as_cortisol_data(data)
  fits <- lapply(alpha_grid, function(a) nmea_fit(data, alpha = a, control = control))
  status <- vapply(fits, `[[`, character(1), "status")
  rss <- vapply(fits, function(f) if (f$status == "ok") f$rss else Inf, numeric(1))
  # Compact diagnostics for every candidate, kept even when keep_fits = FALSE.
  profile <- data.frame(
    alpha = alpha_grid, status = status, rss = rss,
    fim_ok = vapply(fits, function(f) if (f$status == "ok") f$fim_ok else NA, logical(1)),
    n_warnings = vapply(fits, function(f) length(f$warnings), integer(1)),
    warnings = vapply(fits, function(f) paste(f$warnings, collapse = " | "), character(1)),
    error = vapply(fits, function(f) if (is.null(f$error)) NA_character_ else f$error, character(1)),
    elapsed = vapply(fits, function(f) f$elapsed, numeric(1)),
    stringsAsFactors = FALSE
  )
  if (!any(is.finite(rss))) {
    reasons <- unique(stats::na.omit(profile$error))
    return(structure(list(status = "failure",
                          error = paste("All alpha fits failed:", paste(reasons, collapse = " | ")),
                          profile = profile, fits = if (keep_fits) fits),
                     class = "nmea_alpha_fit"))
  }
  best <- which.min(rss)
  structure(list(status = "ok", alpha = alpha_grid[best], best = fits[[best]],
                 profile = profile, fits = if (keep_fits) fits),
            class = "nmea_alpha_fit")
}

#' @export
print.nmea_alpha_fit <- function(x, ...) {
  cat(sprintf("<nmea_alpha_fit> %s; selected alpha = %s\n", x$status,
              if (x$status == "ok") format(x$alpha) else "none"))
  print(x$profile[, c("alpha", "status", "rss", "fim_ok", "n_warnings")], row.names = FALSE)
  invisible(x)
}
