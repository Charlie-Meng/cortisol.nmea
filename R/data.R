#' Create a validated cortisol data set
#'
#' Standardizes a long-format data frame of cortisol measurements for use with
#' the fitting and simulation functions. Rows keep their original order.
#'
#' The time axis is hours since waking; the first sample of a day is taken as
#' the waking time. Concentrations must be in nmol/L (see
#' [convert_cortisol_units()]).
#'
#' @param data A data frame in long format, one row per measurement.
#' @param id,time,value Names of the columns holding the subject identifier,
#'   the time since waking in hours, and the cortisol concentration in nmol/L.
#' @param point_id Optional name of a column holding unique measurement
#'   identifiers. If `NULL`, identifiers are created from the subject and the
#'   row order within subject.
#' @param trimester,sample Optional names of columns holding the pregnancy
#'   trimester and the within-day sample number.
#'
#' @return A data frame of class `cortisol_data` with columns `subject`,
#'   `time`, `y`, `point_id`, and, when supplied, `trimester` and `sample`.
#' @export
#' @examples
#' raw <- data.frame(
#'   pid = rep(c("a", "b", "c"), each = 4),
#'   hours = rep(c(0, 0.5, 6, 12), 3),
#'   conc = c(12, 18, 6, 3, 10, 15, 5, 2, 14, 20, 7, 4)
#' )
#' cortisol_data(raw, id = "pid", time = "hours", value = "conc")
cortisol_data <- function(data, id, time, value, point_id = NULL,
                          trimester = NULL, sample = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  cols <- c(id = id, time = time, value = value)
  optional <- c(point_id = point_id, trimester = trimester, sample = sample)
  missing_cols <- setdiff(c(cols, optional), names(data))
  if (length(missing_cols)) {
    stop("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  out <- data.frame(
    subject = as.character(data[[id]]),
    time = data[[time]],
    y = data[[value]],
    stringsAsFactors = FALSE
  )
  if (!is.numeric(out$time) || !is.numeric(out$y)) {
    stop("Time and value columns must be numeric.", call. = FALSE)
  }
  if (anyNA(out$subject)) stop("Subject identifiers must not be missing.", call. = FALSE)
  if (any(!is.finite(out$time)) || any(!is.finite(out$y))) {
    stop("Time and value columns must be finite (no NA, NaN or Inf).", call. = FALSE)
  }
  if (any(out$time < 0)) stop("Time since waking must be non-negative.", call. = FALSE)
  if (any(out$y < 0)) stop("Cortisol concentrations must be non-negative.", call. = FALSE)

  if (is.null(point_id)) {
    within <- stats::ave(seq_len(nrow(out)), out$subject, FUN = seq_along)
    out$point_id <- sprintf("%s_%03d", out$subject, within)
  } else {
    out$point_id <- as.character(data[[point_id]])
  }
  if (anyNA(out$point_id) || anyDuplicated(out$point_id)) {
    stop("Measurement identifiers must be unique and non-missing.", call. = FALSE)
  }
  if (!is.null(trimester)) out$trimester <- as.character(data[[trimester]])
  if (!is.null(sample)) out$sample <- as.integer(data[[sample]])

  counts <- table(out$subject)
  if (any(counts < 3L)) {
    warning(sprintf("%d subject(s) have fewer than 3 measurements; their curves are weakly identified.",
                    sum(counts < 3L)), call. = FALSE)
  }
  rownames(out) <- NULL
  class(out) <- c("cortisol_data", "data.frame")
  out
}

#' @export
print.cortisol_data <- function(x, n = 6L, ...) {
  counts <- table(x$subject)
  cat(sprintf("<cortisol_data> %d measurements from %d subjects (%s per subject)\n",
              nrow(x), length(counts),
              paste(range(counts), collapse = "-")))
  cat(sprintf("Time since waking: %.2f to %.2f h; concentration: %.2f to %.2f nmol/L\n",
              min(x$time), max(x$time), min(x$y), max(x$y)))
  print(utils::head(as.data.frame(unclass(x), stringsAsFactors = FALSE), n), ...)
  invisible(x)
}

# Accept a cortisol_data object or a data frame with the standard columns.
.as_cortisol_data <- function(x) {
  if (inherits(x, "cortisol_data")) return(x)
  if (is.data.frame(x) && all(c("subject", "time", "y") %in% names(x))) {
    return(cortisol_data(
      x, id = "subject", time = "time", value = "y",
      point_id = if ("point_id" %in% names(x)) "point_id" else NULL,
      trimester = if ("trimester" %in% names(x)) "trimester" else NULL,
      sample = if ("sample" %in% names(x)) "sample" else NULL
    ))
  }
  stop("`data` must be created with cortisol_data() or contain columns subject, time and y.",
       call. = FALSE)
}

# Subset rows while keeping the cortisol_data class.
.subset_cd <- function(x, rows) {
  out <- as.data.frame(unclass(x), stringsAsFactors = FALSE)[rows, , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("cortisol_data", "data.frame")
  out
}

#' Convert cortisol concentrations to nmol/L
#'
#' Uses the molar mass of cortisol, 362.46 g/mol.
#'
#' @param x Numeric concentrations.
#' @param from Unit of `x`: `"ug/dL"`, `"ng/mL"` (equivalently ug/L), or
#'   `"nmol/L"`.
#'
#' @return Numeric concentrations in nmol/L.
#' @export
#' @examples
#' convert_cortisol_units(c(0.5, 1), from = "ug/dL")
convert_cortisol_units <- function(x, from = c("ug/dL", "ng/mL", "nmol/L")) {
  from <- match.arg(from)
  factor <- switch(from,
    "ug/dL" = 1e4 / 362.46,
    "ng/mL" = 1e3 / 362.46,
    "nmol/L" = 1
  )
  x * factor
}

#' Restrict cortisol data to a time window
#'
#' Optional preprocessing. The NMEA features and the evaluation grid use
#' 0 to 18 hours after waking by default.
#'
#' @param x A `cortisol_data` object.
#' @param time_window Length-two numeric vector; measurements outside the
#'   closed interval are removed.
#'
#' @return A `cortisol_data` object.
#' @export
#' @examples
#' raw <- data.frame(id = rep(1:3, each = 3), t = rep(c(0, 8, 20), 3),
#'                   y = c(15, 5, 2, 12, 4, 1, 18, 6, 3))
#' d <- suppressWarnings(cortisol_data(raw, "id", "t", "y"))
#' prep_cortisol(d, time_window = c(0, 18))
prep_cortisol <- function(x, time_window = c(0, 18)) {
  x <- .as_cortisol_data(x)
  if (!is.numeric(time_window) || length(time_window) != 2L || time_window[1] > time_window[2]) {
    stop("`time_window` must be an increasing numeric vector of length 2.", call. = FALSE)
  }
  keep <- x$time >= time_window[1] & x$time <= time_window[2]
  .subset_cd(x, keep)
}
