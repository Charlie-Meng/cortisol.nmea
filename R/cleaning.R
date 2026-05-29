#' Fixed residual outlier rule
#'
#' @param residual Numeric residuals.
#' @param cutoff Absolute residual cutoff.
#'
#' @return Logical vector where `TRUE` indicates flagged observations.
#' @export
detect_outliers_fixed <- function(residual, cutoff = 6) {
  abs(residual) > cutoff
}

#' Time-varying residual outlier rule
#'
#' Estimates a local residual scale `sigma(t)` and flags observations by a
#' standardized residual cutoff.
#'
#' @param time Numeric observation times.
#' @param residual Numeric residuals.
#' @param z_cutoff Absolute standardized residual cutoff.
#' @param span Loess span passed to [estimate_sigma_t()].
#'
#' @return A data frame with `Sigma_t`, `Z`, and `Flag`.
#' @export
detect_outliers_dynamic <- function(time, residual, z_cutoff = 3, span = 0.4) {
  sigma_fun <- estimate_sigma_t(time, residual, span = span)
  sigma_t <- sigma_fun(time)
  z <- residual / pmax(sigma_t, 1e-6)
  data.frame(
    Sigma_t = sigma_t,
    Z = z,
    Flag = abs(z) > z_cutoff
  )
}
