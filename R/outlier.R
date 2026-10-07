#' Outlier rules for the NMEA workflow
#'
#' Rules flag measurements from the initial NMEA fit. A measurement is flagged
#' when its absolute residual exceeds the rule's cutoff.
#'
#' * `outlier_default()` is the package default. For data labelled with a
#'   pregnancy trimester (`"T1"`, `"T2"`, `"T3"`), or when `reference` is
#'   given, it uses `outlier_sd(3.5, "reference")`; otherwise
#'   `outlier_sd(3, "iterative")`, an empirical first-pass screen for clearly
#'   aberrant measurements.
#' * `outlier_fixed(cutoff)` uses a fixed cutoff in nmol/L. `cutoff = 6`
#'   reproduces the original NMEA analysis.
#' * `outlier_sd(k, sd)` uses `cutoff = k * sd`. `sd` can be:
#'   - a number;
#'   - `"reference"`: the noise SD of a simulation reference. By default the
#'     built-in reference ([sim_reference()]) matching the data's `trimester`
#'     column; pass `reference` to use another one. In [sim_run()] the
#'     design's own references are used;
#'   - `"iterative"`: remove measurements with absolute initial residual above
#'     `2.75 x MAD`, refit at the initial alpha, and use that fit's residual
#'     SD (needs one extra fit, whose diagnostics are kept);
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
#' If a rule fails (for example the iterative refit does not converge), the
#' workflow returns a failure for that fit instead of stopping.
#'
#' @param cutoff Positive cutoff in nmol/L.
#' @param k Positive multiplier.
#' @param sd Noise scale; see Details.
#' @param reference For `sd = "reference"`: `NULL` (built-in reference chosen
#'   by the trimester label), an `nmea_reference`, or a noise SD in nmol/L.
#' @param fun Function taking an `nmea_fit` and returning a logical vector.
#'
#' @return An object of class `nmea_outlier_rule`.
#' @export
#' @examples
#' outlier_default()
#' outlier_fixed(6)
#' outlier_sd(k = 3.5, sd = "reference")
#' outlier_sd(k = 2.75, sd = "mad")
outlier_default <- function(reference = NULL) {
  ref <- outlier_sd(3.5, "reference", reference = reference)
  iter <- outlier_sd(3, "iterative")
  use_ref <- function(fit) !is.null(reference) || .has_reference_trimester(fit$data)
  structure(list(
    label = "default: 3.5 x reference SD for T1/T2/T3 data, otherwise 3 x iterative SD",
    spec = list(type = "default", reference = ref$spec$reference),
    cutoff_fun = function(fit) {
      if (use_ref(fit)) ref$cutoff_fun(fit) else iter$cutoff_fun(fit)
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
outlier_sd <- function(k = 3.5, sd, reference = NULL) {
  .check_scalar_number(k, "k", lower = 0)
  if (missing(sd)) stop("`sd` must be supplied.", call. = FALSE)
  ref_sd <- .reference_value(reference)
  sd_fun <- if (is.function(sd)) {
    sd
  } else if (is.character(sd) && length(sd) == 1L) {
    switch(sd,
      reference = function(fit) if (is.null(ref_sd)) .reference_sd(fit) else ref_sd,
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
    # Functions are kept as they are and described when the cache key is
    # computed, so later changes to the variables they use are seen.
    spec = list(type = "sd", k = k, sd = sd, reference = ref_sd),
    cutoff_fun = function(fit) {
      value <- sd_fun(fit)
      .check_scalar_number(as.numeric(value), "sd", lower = 0)
      structure(k * as.numeric(value), diagnostics = attr(value, "diagnostics"))
    }
  ), class = "nmea_outlier_rule")
}

#' @rdname outlier_default
#' @export
outlier_custom <- function(fun) {
  if (!is.function(fun)) stop("`fun` must be a function.", call. = FALSE)
  structure(list(label = "custom rule", spec = list(type = "custom", fun = fun),
                 flag_fun = fun), class = "nmea_outlier_rule")
}

#' @export
print.nmea_outlier_rule <- function(x, ...) {
  cat("<nmea_outlier_rule>", x$label, "\n")
  invisible(x)
}

.reference_value <- function(reference) {
  if (is.null(reference)) return(NULL)
  value <- if (inherits(reference, "nmea_reference")) reference$noise_sd else reference
  .check_scalar_number(value, "reference noise SD", lower = 0)
  value
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
  diagnostics <- list(stage = "iterative_sd", status = refit$status, error = refit$error,
                      fim_ok = if (refit$status == "ok") refit$fim_ok else NA,
                      warnings = refit$warnings, elapsed = refit$elapsed, n_kept = sum(keep))
  if (refit$status != "ok") {
    stop("The iterative SD refit failed: ", refit$error, call. = FALSE)
  }
  structure(refit$sigma, diagnostics = diagnostics)
}

# Rules whose "reference" SD was not fixed when they were built take the
# reference of the simulation design, so fitting and tune_threshold() agree.
.with_reference <- function(rule, reference) {
  if (is.null(rule) || is.numeric(rule) || !is.null(rule$spec$reference)) return(rule)
  if (identical(rule$spec$type, "default")) return(outlier_default(reference = reference))
  if (identical(rule$spec$type, "sd") && identical(rule$spec$sd, "reference")) {
    return(outlier_sd(rule$spec$k, "reference", reference = reference))
  }
  rule
}

# Plain-data description of a function for cache keys: its source and the
# current values of every variable it uses that is not in a package, found
# in its enclosing environments (including the global environment).
.function_signature <- function(fun, depth = 3L) {
  if (is.primitive(fun)) return(list(primitive = deparse(fun)))
  env <- environment(fun)
  used <- unique(unlist(codetools::findGlobals(fun, merge = FALSE)))
  local <- if (!is.null(env) && !identical(env, globalenv()) && !isNamespace(env)) {
    ls(env, all.names = TRUE)
  } else {
    character()
  }
  captured <- list()
  if (!is.null(env) && !isNamespace(env)) {
    for (v in sort(unique(c(used, local)))) {
      if (!exists(v, envir = env)) next
      obj <- get(v, envir = env)
      if (is.function(obj)) {
        fenv <- environment(obj)
        if (is.primitive(obj) || is.null(fenv) || isNamespace(fenv) || identical(fenv, baseenv())) next
        captured[[v]] <- if (depth > 0) .function_signature(obj, depth - 1L) else deparse(obj)
      } else if (!is.environment(obj)) {
        captured[[v]] <- obj
      }
    }
  }
  list(body = deparse(fun), captured = captured)
}

# Replace functions inside a spec by their signatures.
.spec_plain <- function(x) {
  if (is.function(x)) return(.function_signature(x))
  if (is.list(x)) return(lapply(x, .spec_plain))
  x
}

# Apply a rule to an initial fit; returns list(flag, cutoff, diagnostics).
.apply_outlier_rule <- function(rule, fit) {
  if (is.numeric(rule)) rule <- outlier_fixed(rule)
  if (!inherits(rule, "nmea_outlier_rule")) {
    stop("`outlier` must be an nmea_outlier_rule, a number, or NULL.", call. = FALSE)
  }
  diagnostics <- NULL
  if (!is.null(rule$flag_fun)) {
    flag <- rule$flag_fun(fit)
    cutoff <- NA_real_
  } else {
    cutoff <- rule$cutoff_fun(fit)
    diagnostics <- attr(cutoff, "diagnostics")
    cutoff <- as.numeric(cutoff)
    flag <- abs(fit$residual) > cutoff
  }
  if (!is.logical(flag) || length(flag) != nrow(fit$data) || anyNA(flag)) {
    stop("The outlier rule must return one non-missing logical value per measurement.",
         call. = FALSE)
  }
  list(flag = flag, cutoff = cutoff, diagnostics = diagnostics)
}
