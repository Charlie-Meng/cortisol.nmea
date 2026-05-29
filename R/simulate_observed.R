#' Estimate a time-varying residual standard deviation curve
#'
#' Estimates `sigma(t)` with a loess fit to squared residuals. The function is
#' generic and can be used with synthetic residuals or approved public data.
#'
#' @param time Numeric observation times.
#' @param residual Numeric residuals.
#' @param span Loess span.
#' @param min_var Lower bound for variance estimates.
#' @param fallback_sd Fallback standard deviation when loess cannot be fit.
#'
#' @return A function that maps numeric times to estimated standard deviations.
#' @export
estimate_sigma_t <- function(time,
                             residual,
                             span = 0.5,
                             min_var = 1e-4,
                             fallback_sd = stats::sd(residual, na.rm = TRUE)) {
  dat <- data.frame(Time = as.numeric(time), ResSq = as.numeric(residual)^2)
  dat <- dat[is.finite(dat$Time) & is.finite(dat$ResSq), , drop = FALSE]

  if (!is.finite(fallback_sd) || fallback_sd <= 0) fallback_sd <- 1

  if (nrow(dat) < 10 || length(unique(dat$Time)) < 4) {
    return(function(tt) rep(fallback_sd, length(tt)))
  }

  fit <- tryCatch(
    stats::loess(
      ResSq ~ Time,
      data = dat,
      span = span,
      control = stats::loess.control(surface = "direct")
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(function(tt) rep(fallback_sd, length(tt)))
  }

  function(tt) {
    var_hat <- as.numeric(stats::predict(fit, newdata = data.frame(Time = as.numeric(tt))))
    var_hat[!is.finite(var_hat)] <- fallback_sd^2
    sqrt(pmax(var_hat, min_var))
  }
}

#' Simulate observed cortisol measurements from latent NMEA parameters
#'
#' Converts synthetic latent subject-level NMEA curves into observed repeated
#' measurements using sampled time templates and observation noise.
#'
#' @param psi Parameter matrix with columns `mu`, `s`, `c1`, `c0`, `alpha`.
#' @param time_templates Either a list of numeric time vectors or a data frame
#'   containing template identifiers and time values.
#' @param sigma Numeric scalar, numeric function `sigma(time)`, or numeric vector
#'   recycled within each subject.
#' @param seed Optional random seed.
#' @param template_id Column name for template IDs when `time_templates` is a
#'   data frame.
#' @param time Column name for times when `time_templates` is a data frame.
#' @param min_value Lower bound applied to noisy cortisol observations.
#'
#' @return A data frame with synthetic observed data.
#' @export
simulate_observed_data <- function(psi,
                                   time_templates = list(c(0, 0.5, 1, 3, 6, 9, 12, 15)),
                                   sigma = 1,
                                   seed = NULL,
                                   template_id = "Template",
                                   time = "Time",
                                   min_value = 0.01) {
  if (!is.null(seed)) set.seed(seed)
  psi <- .as_psi_matrix(psi, require_alpha = TRUE)

  templates <- if (is.data.frame(time_templates)) {
    if (!all(c(template_id, time) %in% names(time_templates))) {
      stop("`time_templates` must contain `template_id` and `time` columns.", call. = FALSE)
    }
    split(time_templates[[time]], time_templates[[template_id]])
  } else {
    time_templates
  }
  templates <- lapply(templates, function(x) sort(as.numeric(x[is.finite(x)])))
  templates <- templates[lengths(templates) > 0]
  if (length(templates) == 0) {
    stop("At least one non-empty time template is required.", call. = FALSE)
  }

  rows <- vector("list", nrow(psi))
  subject_ids <- rownames(psi)
  if (is.null(subject_ids)) subject_ids <- paste0("SIM_", sprintf("%03d", seq_len(nrow(psi))))

  for (i in seq_len(nrow(psi))) {
    tt <- sample(templates, size = 1)[[1]]
    truth <- asym_sl_model(tt, psi[i, , drop = FALSE])
    sigma_t <- if (is.function(sigma)) sigma(tt) else rep(sigma, length.out = length(tt))
    obs <- pmax(truth + stats::rnorm(length(tt), mean = 0, sd = sigma_t), min_value)

    rows[[i]] <- data.frame(
      Subject = subject_ids[i],
      Time = tt,
      Truth = truth,
      Sigma_t = sigma_t,
      Conc_Obs = obs,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
