#' Simulate subject-level NMEA parameters
#'
#' Draws synthetic subject-level parameters from a multivariate normal
#' distribution with simple biological constraints. This function is intended
#' for public simulation examples; it should not be fed private fitted
#' subject-level estimates unless those estimates are approved for sharing.
#'
#' @param n Number of synthetic subjects.
#' @param theta Mean vector for `mu`, `s`, `c1`, and `c0`.
#' @param covariance Variance-covariance matrix for `mu`, `s`, `c1`, and `c0`.
#' @param alpha Fixed alpha value or length-`n` vector.
#' @param seed Optional random seed.
#' @param min_s Minimum accepted `s`.
#' @param min_c1 Minimum accepted `c1`.
#' @param min_c0 Minimum accepted `c0`.
#' @param max_attempts Maximum rejection-sampling attempts.
#'
#' @return A parameter matrix with columns `mu`, `s`, `c1`, `c0`, `alpha`.
#' @export
simulate_nmea_parameters <- function(n,
                                     theta = c(mu = 1.0, s = 1.2, c1 = 35, c0 = 2),
                                     covariance = diag(c(0.25, 0.12, 5, 0.5)^2),
                                     alpha = 1.15,
                                     seed = NULL,
                                     min_s = 0.05,
                                     min_c1 = 0.05,
                                     min_c0 = 0,
                                     max_attempts = 200) {
  if (!is.null(seed)) set.seed(seed)
  theta <- theta[c("mu", "s", "c1", "c0")]
  covariance <- as.matrix(covariance)
  if (is.null(colnames(covariance)) || is.null(rownames(covariance))) {
    if (!identical(dim(covariance), c(4L, 4L))) {
      stop("`covariance` must be a 4 x 4 matrix for mu, s, c1, and c0.", call. = FALSE)
    }
    colnames(covariance) <- c("mu", "s", "c1", "c0")
    rownames(covariance) <- c("mu", "s", "c1", "c0")
  }
  covariance <- covariance[c("mu", "s", "c1", "c0"), c("mu", "s", "c1", "c0")]
  covariance <- .make_positive_definite(covariance)

  out <- matrix(NA_real_, nrow = n, ncol = 4)
  colnames(out) <- c("mu", "s", "c1", "c0")
  kept <- 0L

  for (attempt in seq_len(max_attempts)) {
    if (kept >= n) break
    batch_n <- max(n - kept, n)
    z <- matrix(stats::rnorm(batch_n * 4), ncol = 4)
    draw <- sweep(z %*% chol(covariance), 2, theta, "+")
    colnames(draw) <- colnames(out)
    ok <- is.finite(rowSums(draw)) &
      draw[, "s"] >= min_s &
      draw[, "c1"] >= min_c1 &
      draw[, "c0"] >= min_c0
    draw <- draw[ok, , drop = FALSE]
    if (nrow(draw) == 0) next

    take <- min(nrow(draw), n - kept)
    out[(kept + 1):(kept + take), ] <- draw[seq_len(take), , drop = FALSE]
    kept <- kept + take
  }

  if (kept < n) {
    stop("Could not draw enough valid NMEA parameter sets.", call. = FALSE)
  }

  alpha <- rep(alpha, length.out = n)
  out <- cbind(out, alpha = alpha)
  rownames(out) <- paste0("SIM_", sprintf("%03d", seq_len(n)))
  out
}

#' Simulate latent NMEA cortisol curves
#'
#' @param psi Parameter matrix with columns `mu`, `s`, `c1`, `c0`, `alpha`.
#' @param time_grid Numeric time grid.
#'
#' @return A list containing parameters, time grid, and curve matrix.
#' @export
simulate_latent_curves <- function(psi, time_grid = seq(0, 18, by = 0.1)) {
  psi <- .as_psi_matrix(psi, require_alpha = TRUE)
  list(
    psi = psi,
    time_grid = time_grid,
    curve_mat = model_eval(time_grid, psi)
  )
}
