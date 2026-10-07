#' GAMM comparator in the style of Sánchez et al. (2012)
#'
#' Fits a generalized additive mixed model to `log1p(cortisol)` on the
#' `sqrt(time / 24)` scale with a penalized regression spline for the
#' population curve, and subject-level random intercepts and slopes plus a
#' subject-level random spline deviation, following the functional mixed model
#' with penalized splines described by Sánchez et al. (2012). Predictions are
#' back-transformed with a conditional log-normal mean correction.
#'
#' This is an adaptation for comparison purposes, **not the original authors'
#' code**. The random spline basis used for prediction at new times is
#' reconstructed from the mgcv model matrix and verified against the fitted
#' values at the observed times.
#'
#' A numerical failure is returned as `status = "failure"` with the error
#' message; no fallback model is substituted.
#'
#' @param data A [cortisol_data()] object.
#' @param times Times since waking (hours) at which to predict curves.
#' @param k Basis dimension of the spline.
#'
#' @return An object of class `gamm_sanchez_fit` with elements `status`,
#'   `curves` (subjects x `times` matrix), `fitted` (predictions at the
#'   observed measurements, in input row order), `times`, `basis_check`,
#'   `warnings`, `error` and `elapsed`.
#' @references Sánchez BN, Wu M, Raghunathan TE, Diez-Roux AV (2012). Modeling
#'   the salivary cortisol profile in population research: the Multi-Ethnic
#'   Study of Atherosclerosis. *American Journal of Epidemiology*, 176(10),
#'   918-928. \doi{10.1093/aje/kws182}
#'
#'   Wood SN (2017). *Generalized Additive Models: An Introduction with R*,
#'   2nd edition. Chapman and Hall/CRC.
#' @export
fit_gamm_sanchez <- function(data, times = seq(0, 18, by = 0.1), k = 5) {
  data <- .as_cortisol_data(data)
  ids <- sort(unique(data$subject))
  started <- proc.time()[["elapsed"]]
  res <- .collect_conditions(.fit_gamm_sanchez(data, ids, times, k))
  elapsed <- proc.time()[["elapsed"]] - started
  if (!is.null(res$error)) {
    return(structure(list(
      status = "failure", error = res$error, times = times,
      curves = matrix(NA_real_, length(ids), length(times), dimnames = list(ids, NULL)),
      fitted = rep(NA_real_, nrow(data)), warnings = res$warnings, elapsed = elapsed
    ), class = "gamm_sanchez_fit"))
  }
  out <- res$value
  out$status <- "ok"
  out$error <- NULL
  out$warnings <- res$warnings
  out$elapsed <- elapsed
  structure(out, class = "gamm_sanchez_fit")
}

.fit_gamm_sanchez <- function(data, ids, times, k) {
  # nlme labels nested random effects "outer/inner", so user identifiers (which
  # may contain "/" or other separators) are replaced by safe internal labels in
  # the same order; the factor coding, and hence the fit, is unchanged.
  internal <- sprintf("s%07d", seq_along(ids))
  dat <- data.frame(subject = factor(internal[match(data$subject, ids)], levels = internal),
                    x = sqrt(data$time / 24), z = log1p(data$y))
  first <- mgcv::gamm(z ~ s(x, k = k), data = dat, method = "REML")
  dat$Xr <- first$lme$data$Xr
  basis <- stats::predict(first$gam, type = "lpmatrix")
  rotation <- qr.solve(basis, as.matrix(dat$Xr))
  if (max(abs(basis %*% rotation - as.matrix(dat$Xr))) >= 1e-7) {
    stop("Random spline basis could not be reconstructed.", call. = FALSE)
  }
  fit <- mgcv::gamm(z ~ s(x, k = k), data = dat,
                    random = list(subject = nlme::pdSymm(~x), subject = nlme::pdIdent(~Xr - 1)),
                    method = "REML")
  re <- fit$lme$coefficients$random
  linear <- re[[2]]
  curved <- re[[3]]
  to_id <- function(labels) ids[match(sub("^.*/", "", labels), internal)]
  rownames(linear) <- to_id(rownames(linear))
  rownames(curved) <- to_id(rownames(curved))
  if (anyNA(rownames(linear)) || anyNA(rownames(curved)) ||
      !setequal(rownames(linear), ids) || !setequal(rownames(curved), ids)) {
    stop("Subject random effects could not be matched.", call. = FALSE)
  }
  sigma <- fit$lme$sigma
  at <- function(subject, time) {
    nd <- data.frame(x = sqrt(time / 24))
    xnew <- stats::predict(first$gam, nd, type = "lpmatrix") %*% rotation
    eta <- as.numeric(stats::predict(fit$gam, nd)) +
      rowSums(linear[subject, , drop = FALSE] * cbind(1, nd$x)) +
      rowSums(curved[subject, , drop = FALSE] * xnew)
    expm1(eta + sigma^2 / 2)
  }
  rows <- at(data$subject, data$time)
  expected <- as.numeric(stats::fitted(fit$lme))
  basis_check <- max(abs(expected - (log1p(rows) - sigma^2 / 2)))
  if (basis_check >= 1e-6) stop("Reconstructed predictions do not match the fitted values.", call. = FALSE)
  curves <- matrix(at(rep(ids, each = length(times)), rep(times, length(ids))),
                   length(ids), byrow = TRUE, dimnames = list(ids, NULL))
  list(curves = curves, fitted = rows, times = times, basis_check = basis_check)
}

#' @export
print.gamm_sanchez_fit <- function(x, ...) {
  if (x$status != "ok") {
    cat("<gamm_sanchez_fit> FAILED:", x$error, "\n")
  } else {
    cat(sprintf("<gamm_sanchez_fit> %d subjects; basis check %.2e\n", nrow(x$curves), x$basis_check))
  }
  invisible(x)
}
