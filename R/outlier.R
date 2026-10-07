#' Outlier rules for the NMEA workflow
#'
#' Rules flag measurements from the initial NMEA fit. A measurement is flagged
#' when its absolute residual exceeds the rule's cutoff.
#'
#' * `outlier_default()` is the package default. For data labelled with a
#'   pregnancy trimester (`"T1"`, `"T2"`, `"T3"`) it uses
#'   `outlier_sd(3.5, "reference")`; otherwise `outlier_sd(3, "iterative")`,
#'   an empirical first-pass screen for clearly aberrant measurements.
#' * `outlier_fixed(cutoff)` uses a fixed cutoff in nmol/L. `cutoff = 6`
#'   reproduces the original NMEA analysis.
#' * `outlier_sd(k, sd)` uses `cutoff = k * sd`. `sd` can be:
#'   - a number;
#'   - `"reference"`: the noise SD of the pregnancy-trimester reference
#'     ([sim_reference()]) matching the data's `trimester` column;
#'   - `"iterative"`: remove measurements with absolute initial residual above
#'     `2.75 x MAD`, refit at the initial alpha, and use that fit's residual
#'     SD (needs one extra fit);
#'   - `"model"`: the residual SD of the initial fit (inflated by
#'     contamination, so the cutoff rises as contamination increases);
#'   - `"mad"`: the scaled median absolute deviation of the initial residuals;
#'   - a function of the initial `nmea_fit` returning a number.
#' * `outlier_custom(fun)` accepts any function of the initial `nmea_fit`
#'   returning one logical value per measurement (in the order of `fit$data`).
#'
#' The multipliers 3.5 (reference) and 3 (iterative) were calibrated so that
#' the cutoff for second-trimester data is close to the original 6 nmol/L, and
#' were compared with other choices in a development simulation.
#'
#' @param cutoff Positive cutoff in nmol/L.
#' @param k Positive multiplier.
#' @param sd Noise scale; see Details.
#' @param fun Function taking an `nmea_fit` and returning a logical vector.
#'
#' @return An object of class `nmea_outlier_rule`.
#' @export
#' @examples
#' outlier_default()
#' outlier_fixed(6)
#' outlier_sd(k = 3.5, sd = "reference")
#' outlier_sd(k = 2.75, sd = "mad")
outlier_default <- function() {
  ref <- outlier_sd(3.5, "reference")
  iter <- outlier_sd(3, "iterative")
  structure(list(
    label = "default: 3.5 x reference SD for T1/T2/T3 data, otherwise 3 x iterative SD",
    spec = list(type = "default", reference = ref$spec, iterative = iter$spec),
    cutoff_fun = function(fit) {
      if (.has_reference_trimester(fit$data)) ref$cutoff_fun(fit) else iter$cutoff_fun(fit)
    }
  ), class = "nmea_outlier_rule")
}

#' @rdname outlier_default
#' @export
outlier_fixed <- function(cutoff = 6) {
  .check_scalar_number(cutoff, "cutoff", lower = 0)
  structure(list(
    label = sprintf("fixed |residual| > %s", format(cutoff)),
    spec = list(type = "fixed", cutoff = cutoff),
    cutoff_fun = function(fit) cutoff
  ), class = "nmea_outlier_rule")
}

#' @rdname outlier_default
#' @export
outlier_sd <- function(k = 3.5, sd) {
  .check_scalar_number(k, "k", lower = 0)
  if (missing(sd)) stop("`sd` must be supplied.", call. = FALSE)
  sd_fun <- if (is.function(sd)) {
    sd
  } else if (is.character(sd) && length(sd) == 1L) {
    switch(sd,
      reference = .reference_sd,
      iterative = .iterative_sd,
      model = function(fit) fit$sigma,
      mad = function(fit) stats::mad(fit$residual),
      stop("`sd` must be a number, \"reference\", \"iterative\", \"model\", \"mad\" or a function.",
           call. = FALSE)
    )
  } else {
    .check_scalar_number(sd, "sd", lower = 0)
    function(fit) sd
  }
  sd_label <- if (is.function(sd)) "function" else format(sd)
  structure(list(
    label = sprintf("|residual| > %s x sd (sd = %s)", format(k), sd_label),
    spec = list(type = "sd", k = k, sd = if (is.function(sd)) .function_signature(sd) else sd),
    cutoff_fun = function(fit) {
      value <- sd_fun(fit)
      .check_scalar_number(value, "sd", lower = 0)
      k * value
    }
  ), class = "nmea_outlier_rule")
}

#' @rdname outlier_default
#' @export
outlier_custom <- function(fun) {
  if (!is.function(fun)) stop("`fun` must be a function.", call. = FALSE)
  structure(list(label = "custom rule", spec = list(type = "custom", fun = .function_signature(fun)),
                 flag_fun = fun), class = "nmea_outlier_rule")
}

#' @export
print.nmea_outlier_rule <- function(x, ...) {
  cat("<nmea_outlier_rule>", x$label, "\n")
  invisible(x)
}

.has_reference_trimester <- function(data) {
  "trimester" %in% names(data) && length(unique(data$trimester)) == 1L &&
    unique(data$trimester) %in% c("T1", "T2", "T3")
}

.reference_sd <- function(fit) {
  if (!.has_reference_trimester(fit$data)) {
    stop("sd = \"reference\" needs data with a single trimester label T1, T2 or T3.", call. = FALSE)
  }
  sim_reference(unique(fit$data$trimester))$noise_sd
}

.iterative_sd <- function(fit, k_screen = 2.75) {
  keep <- abs(fit$residual) <= k_screen * stats::mad(fit$residual)
  control <- if (is.null(fit$control)) nmea_control() else fit$control
  refit <- nmea_fit(.subset_cd(fit$data, keep), alpha = fit$alpha, control = control)
  if (refit$status != "ok") stop("The iterative SD refit failed: ", refit$error, call. = FALSE)
  refit$sigma
}

# Plain-data description of a function, including the values it captured, so
# that two closures with the same body but different captured values differ.
.function_signature <- function(fun) {
  env <- environment(fun)
  captured <- if (is.null(env) || identical(env, globalenv()) || isNamespace(env)) {
    NULL
  } else {
    vals <- mget(ls(env, all.names = TRUE), envir = env)
    lapply(vals, function(v) if (is.function(v)) deparse(v) else v)
  }
  list(body = deparse(fun), captured = captured)
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
