#' Asymmetric time warping of the NMEA model
#'
#' Before the peak time `mu`, time is compressed (alpha > 1) or stretched
#' (alpha < 1) towards `mu`; after `mu` it is unchanged:
#' \deqn{t^* = \alpha t + (1 - \alpha)\mu \quad (t < \mu), \qquad t^* = t \quad (t \ge \mu).}
#'
#' @param t Numeric times since waking (hours), or warped times for
#'   `nmea_unwarp_time()`.
#' @param mu Peak location; scalar or the same length as `t`.
#' @param alpha Asymmetry parameter (positive); scalar or the same length as `t`.
#'
#' @return Numeric vector of warped (or unwarped) times.
#' @export
#' @examples
#' tt <- c(0, 0.25, 0.5, 2)
#' w <- nmea_warp_time(tt, mu = 0.5, alpha = 1.5)
#' nmea_unwarp_time(w, mu = 0.5, alpha = 1.5)
nmea_warp_time <- function(t, mu, alpha) {
  tt <- alpha * t + (1 - alpha) * mu
  keep <- which(t >= mu)
  tt[keep] <- t[keep]
  tt
}

#' @rdname nmea_warp_time
#' @export
nmea_unwarp_time <- function(t, mu, alpha) {
  tt <- (t - (1 - alpha) * mu) / alpha
  keep <- which(t >= mu)
  tt[keep] <- t[keep]
  tt
}

#' Evaluate NMEA cortisol curves
#'
#' The asymmetric scaled-logistic curve
#' \deqn{f(t) = c_0 + \frac{c_1}{s} g\left(\frac{t^* - \mu}{s}\right),}
#' where \eqn{g} is the standard logistic density and \eqn{t^*} the warped time
#' from [nmea_warp_time()]. `mu` is the peak time, `s` the width, `c1` the
#' amplitude (a positive `c1` gives a morning peak), `c0` the baseline and
#' `alpha` the pre-peak asymmetry.
#'
#' The logistic density is evaluated in a numerically stable form, so very
#' large `|t - mu| / s` gives the baseline instead of overflowing.
#'
#' @param t Numeric times since waking (hours).
#' @param psi Named numeric vector, or matrix with one row per subject and
#'   columns `mu`, `s`, `c1`, `c0` and `alpha`.
#'
#' @return A numeric vector (for a parameter vector) or a matrix with one row
#'   per subject and one column per time.
#' @export
#' @examples
#' psi <- c(mu = 0.6, s = 0.75, c1 = 30, c0 = 2.5, alpha = 1.5)
#' nmea_curve(c(0, 0.5, 1, 6, 12), psi)
nmea_curve <- function(t, psi) {
  single <- is.null(dim(psi))
  p <- .as_psi_matrix(psi)
  out <- matrix(NA_real_, nrow(p), length(t), dimnames = list(rownames(p), NULL))
  for (i in seq_len(nrow(p))) {
    out[i, ] <- .nmea_curve_one(t, p[i, ])
  }
  if (single) out[1, ] else out
}

.nmea_curve_one <- function(t, theta) {
  s <- theta[["s"]]
  if (!is.finite(s) || abs(s) < 1e-12) return(rep(NA_real_, length(t)))
  tt <- nmea_warp_time(t, theta[["mu"]], theta[["alpha"]])
  theta[["c0"]] + theta[["c1"]] / s * stats::dlogis((tt - theta[["mu"]]) / s)
}

# Structural model handed to saemix. This must stay arithmetically identical to
# the original NMEA_Cortisol.R `modelFun` (same operations in the same order),
# because SAEM iterations are sensitive to last-bit floating-point differences
# and published results must be reproducible exactly.
.nmea_saemix_model <- function(psi, id, x) {
  tp <- x[, 1]
  mu <- psi[id, 1]
  s <- psi[id, 2]
  c1 <- psi[id, 3]
  c0 <- psi[id, 4]
  alpha <- psi[id, 5]
  tp.new <- alpha * tp + (1 - alpha) * mu
  nochange <- which(tp >= mu)
  tp.new[nochange] <- tp[nochange]
  c1 * exp(-(tp.new - mu) / s) / (s * (1 + exp(-(tp.new - mu) / s))^2) + c0
}
