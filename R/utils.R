# Internal helpers shared across the package.

# Evaluate `expr` without leaving a changed global RNG state behind.
# saemix seeds the global RNG internally; users' streams should not be altered.
.with_preserved_rng <- function(expr) {
  env <- globalenv()
  had_seed <- exists(".Random.seed", envir = env, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = env, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  }, add = TRUE)
  expr
}

# Evaluate `expr`, discarding console output when `quiet` is TRUE.
.maybe_quiet <- function(expr, quiet = TRUE) {
  if (!isTRUE(quiet)) return(expr)
  out <- NULL
  # saemix also reports internal try() failures on stderr; capture both streams.
  utils::capture.output(utils::capture.output(out <- suppressMessages(expr), type = "message"))
  out
}

# Evaluate `expr`, collecting (and muffling) warnings and catching errors.
# Returns list(value, warnings, error).
.collect_conditions <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) structure(conditionMessage(e), class = ".cn_error")
  )
  if (inherits(value, ".cn_error")) {
    return(list(value = NULL, warnings = unique(warnings), error = unclass(value)))
  }
  list(value = value, warnings = unique(warnings), error = NULL)
}

.check_scalar_number <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper) {
    stop(sprintf("`%s` must be a single finite number in [%s, %s].", name, lower, upper),
         call. = FALSE)
  }
  invisible(x)
}

.psi_names <- c("mu", "s", "c1", "c0", "alpha")

# Coerce parameters to a numeric matrix with columns mu, s, c1, c0, alpha.
.as_psi_matrix <- function(psi) {
  if (is.null(dim(psi))) {
    if (is.null(names(psi)) || !all(.psi_names %in% names(psi))) {
      stop("`psi` must be named with mu, s, c1, c0 and alpha.", call. = FALSE)
    }
    psi <- matrix(psi[.psi_names], nrow = 1, dimnames = list(NULL, .psi_names))
  } else {
    psi <- as.matrix(psi)
    if (is.null(colnames(psi)) || !all(.psi_names %in% colnames(psi))) {
      stop("`psi` must have columns mu, s, c1, c0 and alpha.", call. = FALSE)
    }
    psi <- psi[, .psi_names, drop = FALSE]
  }
  storage.mode(psi) <- "double"
  psi
}
