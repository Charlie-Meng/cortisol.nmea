#' Fit the NMEA model with saemix
#'
#' Fits the NMEA nonlinear mixed-effects model to repeated cortisol
#' measurements. The `alpha` asymmetry parameter is fixed for a given fit;
#' use [fit_alpha_grid()] to search over candidate alpha values.
#'
#' @param data Data frame containing subject, time, and response columns.
#' @param id Name of the subject identifier column.
#' @param time Name of the time column.
#' @param response Name of the cortisol response column.
#' @param alpha Fixed asymmetry parameter.
#' @param psi0 Named initial values for `mu`, `s`, `c1`, and `c0`.
#' @param seed Random seed passed to `saemix`.
#' @param units Optional `saemixData()` units list.
#' @param control Optional list merged into the default `saemix` control list.
#' @param quiet Logical; if `TRUE`, suppresses verbose `saemix` console output.
#' @param ... Additional arguments passed to `saemix::saemixModel()`.
#'
#' @return A list with subject-level parameters, fitted values, residuals, RSS,
#'   FVU, and the raw `saemix` fit.
#' @export
fit_nmea <- function(data,
                     id = "ID",
                     time = "Time",
                     response = "Cortisol",
                     alpha = 1,
                     psi0 = c(mu = 1, s = 1, c1 = 1, c0 = 0),
                     seed = 632545,
                     units = list(x = "Hour", y = "nmol/L"),
                     control = list(),
                     quiet = TRUE,
                     ...) {
  if (!requireNamespace("saemix", quietly = TRUE)) {
    stop("Package `saemix` is required to fit NMEA models.", call. = FALSE)
  }

  data <- as.data.frame(data)
  required <- c(id, time, response)
  if (!all(required %in% names(data))) {
    stop("`data` must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }

  y <- data[[response]]
  subject_id <- as.character(data[[id]])
  model_data <- .quiet_eval(
    saemix::saemixData(
      name.data = data,
      name.group = id,
      name.predictors = time,
      name.response = response,
      units = units
    ),
    quiet = quiet
  )

  psi2 <- c(psi0, alpha = alpha)
  covariance_model <- diag(5)
  covariance_model[5, 5] <- 0

  model <- .quiet_eval(
    saemix::saemixModel(
      model = .saemix_model_fun,
      psi0 = psi2,
      fixed.estim = c(1, 1, 1, 1, 0),
      covariance.model = covariance_model,
      ...
    ),
    quiet = quiet
  )

  default_control <- list(
    map = TRUE,
    fim = TRUE,
    ll.is = FALSE,
    displayProgress = FALSE,
    seed = seed,
    print = FALSE,
    save = FALSE,
    save.graphs = FALSE
  )
  fit <- .quiet_eval(
    saemix::saemix(
      model,
      model_data,
      control = utils::modifyList(default_control, control)
    ),
    quiet = quiet
  )

  psi_est <- saemix::psi(fit)
  rownames(psi_est) <- unique(subject_id)
  ypred <- .predict_observed_rows(
    data = data,
    id = id,
    time = time,
    psi = psi_est
  )
  residual <- y - ypred
  rss <- sum(residual^2, na.rm = TRUE)

  ss_err <- tapply((y - ypred)^2, subject_id, sum, na.rm = TRUE)
  ss_tot <- tapply(y, subject_id, function(yi) sum((yi - mean(yi, na.rm = TRUE))^2, na.rm = TRUE))
  fvu <- ss_err / ss_tot

  structure(
    list(
      psi = psi_est,
      ypred = ypred,
      residual = residual,
      RSS = rss,
      fvu = fvu,
      fit = fit,
      alpha = alpha,
      columns = list(id = id, time = time, response = response)
    ),
    class = "nmea_fit"
  )
}

#' Fit NMEA over an alpha grid
#'
#' @inheritParams fit_nmea
#' @param alpha_grid Numeric vector of candidate alpha values.
#'
#' @return A list containing all successful fits, RSS values, the selected
#'   `alpha_star`, and the final fit.
#' @export
fit_alpha_grid <- function(data,
                           alpha_grid = seq(0.8, 1.5, by = 0.1),
                           id = "ID",
                           time = "Time",
                           response = "Cortisol",
                           psi0 = c(mu = 1, s = 1, c1 = 1, c0 = 0),
                           seed = 632545,
                           quiet = TRUE,
                           ...) {
  fits <- vector("list", length(alpha_grid))
  names(fits) <- paste0("a=", alpha_grid)
  rss <- rep(NA_real_, length(alpha_grid))
  names(rss) <- names(fits)

  for (i in seq_along(alpha_grid)) {
    fits[[i]] <- tryCatch(
      fit_nmea(
        data = data,
        id = id,
        time = time,
        response = response,
        alpha = alpha_grid[i],
        psi0 = psi0,
        seed = seed,
        quiet = quiet,
        ...
      ),
      error = function(e) e
    )
    if (!inherits(fits[[i]], "error")) {
      rss[i] <- fits[[i]]$RSS
    }
  }

  if (!any(is.finite(rss))) {
    stop("No alpha-grid fit succeeded.", call. = FALSE)
  }

  best <- which.min(rss)
  list(
    fits = fits,
    rss = rss,
    alpha_star = alpha_grid[best],
    final_fit = fits[[best]]
  )
}

#' Predict NMEA curves from fitted or supplied parameters
#'
#' @param object An `nmea_fit` object returned by [fit_nmea()] or a parameter
#'   matrix with columns `mu`, `s`, `c1`, `c0`, and `alpha`.
#' @param time_grid Numeric vector of times since waking.
#'
#' @return A numeric prediction matrix.
#' @export
predict_nmea <- function(object, time_grid) {
  psi <- if (inherits(object, "nmea_fit")) object$psi else object
  model_eval(time_grid, psi)
}

.predict_observed_rows <- function(data, id, time, psi) {
  psi <- .as_psi_matrix(psi, require_alpha = TRUE)
  subject_id <- as.character(data[[id]])
  pred <- rep(NA_real_, nrow(data))
  common <- intersect(unique(subject_id), rownames(psi))
  for (sid in common) {
    idx <- which(subject_id == sid)
    pred[idx] <- asym_sl_model(data[[time]][idx], psi[sid, , drop = FALSE])
  }
  pred
}

.quiet_eval <- function(expr, quiet = TRUE) {
  if (!isTRUE(quiet)) return(expr)
  out <- NULL
  utils::capture.output(out <- expr)
  out
}
