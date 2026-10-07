#' Outlier rules for the NMEA workflow
#'
#' Rules flag measurements from the initial NMEA fit. A measurement is flagged
#' when its absolute residual exceeds the rule's cutoff.
#'
#' * `outlier_fixed(cutoff)` uses a fixed cutoff in nmol/L. `cutoff = 6`
#'   reproduces the original NMEA analysis.
#' * `outlier_sd(k, sd)` uses `cutoff = k * sd`, which adapts the cutoff to the
#'   noise level of the data. `sd` can be a number, `"model"` (the residual
#'   error SD estimated by the initial fit), `"mad"` (the scaled median
#'   absolute deviation of the initial residuals), or a function of the
#'   initial `nmea_fit` returning a number.
#' * `outlier_custom(fun)` accepts any function of the initial `nmea_fit`
#'   returning one logical value per measurement (in the order of `fit$data`).
#'
#' @param cutoff Positive cutoff in nmol/L.
#' @param k Positive multiplier.
#' @param sd Noise scale; see Details.
#' @param fun Function taking an `nmea_fit` and returning a logical vector.
#'
#' @return An object of class `nmea_outlier_rule`.
#' @export
#' @examples
#' outlier_fixed(6)
#' outlier_sd(k = 3.5, sd = 1.77)
#' outlier_sd(k = 3.5, sd = "mad")
outlier_fixed <- function(cutoff = 6) {
  .check_scalar_number(cutoff, "cutoff", lower = 0)
  structure(list(
    label = sprintf("fixed |residual| > %s", format(cutoff)),
    spec = list(type = "fixed", cutoff = cutoff),
    cutoff_fun = function(fit) cutoff
  ), class = "nmea_outlier_rule")
}

#' @rdname outlier_fixed
#' @export
outlier_sd <- function(k = 3.5, sd) {
  .check_scalar_number(k, "k", lower = 0)
  if (missing(sd)) stop("`sd` must be supplied.", call. = FALSE)
  sd_fun <- if (is.function(sd)) {
    sd
  } else if (is.character(sd) && length(sd) == 1L) {
    switch(sd,
      model = function(fit) fit$sigma,
      mad = function(fit) stats::mad(fit$residual),
      stop("`sd` must be a number, \"model\", \"mad\" or a function.", call. = FALSE)
    )
  } else {
    .check_scalar_number(sd, "sd", lower = 0)
    function(fit) sd
  }
  sd_label <- if (is.function(sd)) "function" else format(sd)
  structure(list(
    label = sprintf("|residual| > %s x sd (sd = %s)", format(k), sd_label),
    spec = list(type = "sd", k = k, sd = if (is.function(sd)) deparse(sd) else sd),
    cutoff_fun = function(fit) {
      value <- sd_fun(fit)
      .check_scalar_number(value, "sd", lower = 0)
      k * value
    }
  ), class = "nmea_outlier_rule")
}

#' @rdname outlier_fixed
#' @export
outlier_custom <- function(fun) {
  if (!is.function(fun)) stop("`fun` must be a function.", call. = FALSE)
  structure(list(label = "custom rule", spec = list(type = "custom", fun = deparse(fun)),
                 flag_fun = fun), class = "nmea_outlier_rule")
}

#' @export
print.nmea_outlier_rule <- function(x, ...) {
  cat("<nmea_outlier_rule>", x$label, "\n")
  invisible(x)
}

# Apply a rule to an initial fit; returns list(flag, cutoff).
.apply_outlier_rule <- function(rule, fit) {
  if (is.numeric(rule)) rule <- outlier_fixed(rule)
  if (!inherits(rule, "nmea_outlier_rule")) {
    stop("`outlier` must be an nmea_outlier_rule, a number, or NULL.", call. = FALSE)
  }
  if (!is.null(rule$flag_fun)) {
    flag <- rule$flag_fun(fit)
    cutoff <- NA_real_
  } else {
    cutoff <- rule$cutoff_fun(fit)
    flag <- abs(fit$residual) > cutoff
  }
  if (!is.logical(flag) || length(flag) != nrow(fit$data) || anyNA(flag)) {
    stop("The outlier rule must return one non-missing logical value per measurement.",
         call. = FALSE)
  }
  list(flag = flag, cutoff = cutoff)
}
