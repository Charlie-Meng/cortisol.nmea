#' Methods compared in a simulation study
#'
#' Method specifications for [sim_design()].
#'
#' * `method_nmea(steps)` runs [nmea_pipeline()] with the given steps; use
#'   [nmea_steps()] for the full workflow and [nmea_steps_direct()] for an
#'   uncleaned fit.
#' * `method_gamm(k)` runs [fit_gamm_sanchez()], an adaptation of the
#'   approach of Sánchez et al. (2012), not the original authors' code.
#' * `default_methods()` returns the three methods of the paper: Direct NMEA,
#'   Full NMEA and GAMM.
#' * `methods_threshold_grid(k, sd)` returns one full-workflow method per value
#'   of `k`, for comparing thresholds with refitting.
#'
#' @param steps An [nmea_steps()] object.
#' @param control An [nmea_control()] object.
#' @param k For `method_gamm()`, the spline basis dimension; for
#'   `methods_threshold_grid()`, the multipliers to compare.
#' @param sd Noise scale passed to [outlier_sd()].
#' @param full_steps Steps for the full workflow in `default_methods()`.
#'
#' @return A method specification (class `nmea_sim_method`) or a named list
#'   of them.
#' @export
#' @examples
#' default_methods()
#' methods_threshold_grid(c(3, 3.5, 4), sd = 1.8)
method_nmea <- function(steps = nmea_steps(), control = nmea_control()) {
  if (!inherits(steps, "nmea_steps")) stop("`steps` must come from nmea_steps().", call. = FALSE)
  structure(list(type = "nmea", steps = steps, control = control), class = "nmea_sim_method")
}

#' @rdname method_nmea
#' @export
method_gamm <- function(k = 5) {
  .check_scalar_number(k, "k", lower = 3)
  structure(list(type = "gamm", k = k), class = "nmea_sim_method")
}

#' @rdname method_nmea
#' @export
default_methods <- function(full_steps = nmea_steps(), control = nmea_control()) {
  list(Direct = method_nmea(nmea_steps_direct(), control),
       Full = method_nmea(full_steps, control),
       GAMM = method_gamm())
}

#' @rdname method_nmea
#' @export
methods_threshold_grid <- function(k, sd, control = nmea_control()) {
  out <- lapply(k, function(kk) method_nmea(nmea_steps(outlier = outlier_sd(kk, sd)), control))
  names(out) <- paste0("Full_k", format(k))
  out
}

#' @export
print.nmea_sim_method <- function(x, ...) {
  if (x$type == "gamm") {
    cat(sprintf("<nmea_sim_method> GAMM (Sanchez et al. 2012 adaptation), k = %s\n", x$k))
  } else {
    cat("<nmea_sim_method> NMEA workflow, outlier rule:",
        if (is.null(x$steps$outlier)) "none" else x$steps$outlier$label, "\n")
  }
  invisible(x)
}

# Fit one method to one simulated data set; returns a slim, serializable result.
.run_method <- function(method, sim, times) {
  ids <- rownames(sim$truth$psi)
  started <- proc.time()[["elapsed"]]
  if (method$type == "gamm") {
    g <- fit_gamm_sanchez(sim$observed, times = times, k = method$k)
    curves <- g$curves[ids, , drop = FALSE]
    out <- list(status = g$status, error = g$error, curves = curves,
                outcome = stats::setNames(ifelse(apply(curves, 1, function(x) all(is.finite(x))),
                                          "retained", "not_fitted"), ids),
                points = NULL, cutoff = NA_real_, alpha = NA_real_, initial_sigma = NA_real_)
  } else {
    r <- nmea_pipeline(sim$observed, steps = method$steps, control = method$control)
    ini <- r$stages$initial
    out <- list(status = r$status, error = r$error,
                curves = stats::predict(r, times = times, subjects = ids),
                outcome = stats::setNames(r$subjects$outcome, r$subjects$subject)[ids],
                points = r$points[, c("point_id", "initial_residual", "flagged")],
                cutoff = r$cutoff, alpha = r$alpha,
                initial_sigma = if (!is.null(ini) && ini$status == "ok") ini$sigma else NA_real_)
  }
  out$minutes <- (proc.time()[["elapsed"]] - started) / 60
  out
}
