.psi_columns <- c("mu", "s", "c1", "c0", "alpha")

.as_psi_matrix <- function(psi, require_alpha = TRUE) {
  if (is.null(psi)) {
    stop("`psi` must not be NULL.", call. = FALSE)
  }

  if (is.vector(psi) && !is.list(psi)) {
    psi_names <- names(psi)
    psi <- matrix(psi, nrow = 1)
    if (!is.null(psi_names)) {
      colnames(psi) <- psi_names
    }
  }

  psi <- as.matrix(psi)

  needed <- if (require_alpha) .psi_columns else .psi_columns[1:4]
  if (!all(needed %in% colnames(psi))) {
    stop(
      "`psi` must contain columns: ",
      paste(needed, collapse = ", "),
      call. = FALSE
    )
  }

  psi[, needed, drop = FALSE]
}

.make_positive_definite <- function(sigma, eps = 1e-6) {
  sigma <- as.matrix(sigma)
  sigma <- (sigma + t(sigma)) / 2
  eig <- eigen(sigma, symmetric = TRUE)
  values <- pmax(eig$values, eps)
  out <- eig$vectors %*% diag(values, nrow = length(values)) %*% t(eig$vectors)
  out <- (out + t(out)) / 2
  colnames(out) <- colnames(sigma)
  rownames(out) <- rownames(sigma)
  out
}

.trapz <- function(y, x) {
  if (length(y) != length(x)) {
    stop("`y` and `x` must have the same length.", call. = FALSE)
  }
  sum((y[-length(y)] + y[-1]) / 2 * diff(x))
}

.ensure_min_keep <- function(keep, min_keep = 3) {
  keep <- as.logical(keep)
  if (sum(keep, na.rm = TRUE) >= min_keep) return(keep)

  candidates <- which(!keep)
  need <- min(min_keep - sum(keep, na.rm = TRUE), length(candidates))
  if (need > 0) {
    keep[sample(candidates, size = need)] <- TRUE
  }
  keep
}
