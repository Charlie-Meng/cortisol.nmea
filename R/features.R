.f0 <- function(x) exp(-x) / ((1 + exp(-x))^2)
.F0 <- function(x) 1 / (1 + exp(-x))

#' Area under an NMEA cortisol curve
#'
#' Computes the analytic NMEA area under the curve from time 0 to `t_max`.
#'
#' @param psi Matrix or data frame with columns `mu`, `s`, `c1`, `c0`, `alpha`.
#' @param t_max Upper integration limit, in hours since waking.
#'
#' @return Numeric vector of AUC values.
#' @export
auc_nmea <- function(psi, t_max = 18) {
  psi <- .as_psi_matrix(psi, require_alpha = TRUE)
  mu <- psi[, "mu"]
  s <- psi[, "s"]
  c0 <- psi[, "c0"]
  c1 <- psi[, "c1"]
  alpha <- psi[, "alpha"]

  out <- c1 / alpha *
    ((1 - alpha) / 2 + alpha * .F0((t_max - mu) / s) - .F0(-alpha * mu / s)) +
    c0 * t_max
  names(out) <- rownames(psi)
  out
}

#' Estimated morning level from NMEA parameters
#'
#' @param psi Matrix or data frame with columns `mu`, `s`, `c1`, `c0`, `alpha`.
#'
#' @return Numeric vector of estimated morning level values.
#' @export
eml_nmea <- function(psi) {
  psi <- .as_psi_matrix(psi, require_alpha = TRUE)
  mu <- psi[, "mu"]
  s <- psi[, "s"]
  c0 <- psi[, "c0"]
  c1 <- psi[, "c1"]
  alpha <- psi[, "alpha"]

  f0vec <- .f0(-alpha * mu / s)
  idx <- which(mu <= 0)
  f0vec[idx] <- .f0(-mu[idx] / s[idx])
  out <- c1 * f0vec / s + c0
  names(out) <- rownames(psi)
  out
}

#' Derived NMEA cortisol summary features
#'
#' Computes the current package-level feature set used in the cortisol NMEA
#' workflow: AUC, EML, PCL, AR, and DDC.
#'
#' @param psi Matrix or data frame with columns `mu`, `s`, `c1`, `c0`, `alpha`.
#' @param t_max Upper integration limit for AUC.
#'
#' @return A data frame with one row per subject.
#' @export
summary_features <- function(psi, t_max = 18) {
  psi <- .as_psi_matrix(psi, require_alpha = TRUE)
  mu <- psi[, "mu"]
  s <- psi[, "s"]
  c0 <- psi[, "c0"]
  c1 <- psi[, "c1"]

  eml <- eml_nmea(psi)
  pcl <- c1 / (4 * s) + c0
  out <- data.frame(
    Subject = rownames(psi),
    AUC = as.numeric(auc_nmea(psi, t_max = t_max)),
    EML = as.numeric(eml),
    PCL = as.numeric(pcl),
    AR = as.numeric((pcl - eml) * (mu > 0)),
    DDC = as.numeric((-sqrt(3) / 18) * c1 / s^2),
    row.names = NULL
  )
  out
}
