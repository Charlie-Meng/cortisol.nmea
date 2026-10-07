#' Interpretable features of NMEA cortisol curves
#'
#' Computes per-subject curve features from NMEA parameters.
#'
#' * `AUC`: area under the curve from waking to `tmax` hours (nmol*h/L).
#' * `EML`: cortisol level at waking, \eqn{f(0)}.
#' * `PCL`: peak cortisol level \eqn{c_1 / (4s) + c_0}, the curve value at `mu`.
#' * `AR`: awakening rise `PCL - EML` when the peak occurs after waking
#'   (`mu > 0`), otherwise 0.
#' * `DDC`: steepest post-peak decline slope \eqn{-(\sqrt{3}/18) c_1 / s^2}
#'   (nmol/L per hour).
#'
#' The AUC integrates the warped curve piecewise (before and after `mu`), so
#' it is exact for any peak location, including peaks before waking.
#'
#' @param x A parameter matrix (columns `mu`, `s`, `c1`, `c0`, `alpha`), a
#'   named parameter vector, or a fitted object from [nmea_fit()] or
#'   [nmea_pipeline()].
#' @param tmax Upper limit of the AUC integral in hours.
#'
#' @return A matrix with one row per subject and columns `AUC`, `EML`, `PCL`,
#'   `AR`, `DDC`.
#' @export
#' @examples
#' psi <- rbind(a = c(mu = 0.6, s = 0.75, c1 = 30, c0 = 2.5, alpha = 1.5),
#'              b = c(mu = -0.3, s = 0.9, c1 = 25, c0 = 2, alpha = 1.5))
#' nmea_features(psi)
nmea_features <- function(x, tmax = 18) {
  .check_scalar_number(tmax, "tmax", lower = 0)
  psi <- if (inherits(x, c("nmea_fit", "nmea_result"))) stats::coef(x) else x
  p <- .as_psi_matrix(psi)
  mu <- p[, "mu"]; s <- p[, "s"]; c1 <- p[, "c1"]; c0 <- p[, "c0"]; alpha <- p[, "alpha"]

  # Warped part on [0, m] and unwarped part on [m, tmax], with m = mu clipped to [0, tmax].
  m <- pmin(pmax(mu, 0), tmax)
  warped <- c1 / alpha * (stats::plogis(alpha * (m - mu) / s) - stats::plogis(-alpha * mu / s))
  unwarped <- c1 * (stats::plogis((tmax - mu) / s) - stats::plogis((m - mu) / s))
  z0 <- ifelse(mu <= 0, -mu / s, -alpha * mu / s)
  eml <- c1 * stats::dlogis(z0) / s + c0
  pcl <- c1 / (4 * s) + c0
  out <- cbind(
    AUC = warped + unwarped + c0 * tmax,
    EML = eml,
    PCL = pcl,
    AR = (pcl - eml) * (mu > 0),
    DDC = (-sqrt(3) / 18) * c1 / s^2
  )
  rownames(out) <- rownames(p)
  out
}
