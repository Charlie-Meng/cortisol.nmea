#' Symmetric scaled logistic cortisol curve
#'
#' Evaluates the NMEA scaled logistic curve without the asymmetry time
#' transformation.
#'
#' @param time Numeric vector of times since waking.
#' @param psi Named numeric vector or matrix with columns `mu`, `s`, `c1`, `c0`.
#'
#' @return Numeric vector of predicted cortisol values.
#' @export
sl_model <- function(time, psi) {
  psi <- .as_psi_matrix(psi, require_alpha = FALSE)
  if (nrow(psi) != 1) {
    stop("`sl_model()` expects exactly one parameter row.", call. = FALSE)
  }

  mu <- psi[1, "mu"]
  s <- psi[1, "s"]
  c1 <- psi[1, "c1"]
  c0 <- psi[1, "c0"]
  c1 * exp(-(time - mu) / s) / (s * (1 + exp(-(time - mu) / s))^2) + c0
}

#' Asymmetric scaled logistic cortisol curve
#'
#' Evaluates the NMEA curve after applying the pre-peak time transformation
#' controlled by `alpha`.
#'
#' @param time Numeric vector of times since waking.
#' @param psi Named numeric vector or matrix with columns `mu`, `s`, `c1`,
#'   `c0`, and `alpha`.
#'
#' @return Numeric vector of predicted cortisol values.
#' @export
asym_sl_model <- function(time, psi) {
  psi <- .as_psi_matrix(psi, require_alpha = TRUE)
  if (nrow(psi) != 1) {
    stop("`asym_sl_model()` expects exactly one parameter row.", call. = FALSE)
  }

  time_new <- time_forward(time, alpha = psi[1, "alpha"], mu = psi[1, "mu"])
  sl_model(time_new, psi)
}

#' Evaluate NMEA curves on a common time grid
#'
#' @param time_grid Numeric vector of times since waking.
#' @param psi Matrix or data frame with one row per subject and columns `mu`,
#'   `s`, `c1`, `c0`, and `alpha`.
#'
#' @return A numeric matrix with one row per subject and one column per time.
#' @export
model_eval <- function(time_grid, psi) {
  psi <- .as_psi_matrix(psi, require_alpha = TRUE)
  out <- t(vapply(
    seq_len(nrow(psi)),
    function(i) asym_sl_model(time_grid, psi[i, , drop = FALSE]),
    numeric(length(time_grid))
  ))

  rownames(out) <- rownames(psi)
  colnames(out) <- paste0("t=", time_grid)
  out
}

.saemix_model_fun <- function(psi, id, x) {
  time <- x[, 1]
  mu <- psi[id, 1]
  s <- psi[id, 2]
  c1 <- psi[id, 3]
  c0 <- psi[id, 4]
  alpha <- psi[id, 5]
  time_new <- time_forward(time, alpha = alpha, mu = mu)
  c1 * exp(-(time_new - mu) / s) /
    (s * (1 + exp(-(time_new - mu) / s))^2) + c0
}
