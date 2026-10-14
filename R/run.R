#' Design a simulation study
#'
#' Specifies trimesters, cohort size, seeds, contamination scenarios and the
#' methods to compare. Each seed defines one base cohort per trimester; all
#' scenarios of that seed share its subjects, sampling times and clean noise.
#'
#' @param trimesters Trimesters to simulate (`"T1"`, `"T2"`, `"T3"`).
#' @param n Subjects per cohort.
#' @param seeds Either a named list with one integer vector of seeds per
#'   trimester (used as given), or an integer vector. A vector is expanded into
#'   distinct trimester-specific seeds `3 * seed + i` for the `i`-th of
#'   T1/T2/T3, so cohorts of different trimesters are independent, as assumed
#'   by the pooled Monte Carlo intervals of [eval_summary()].
#' @param scenarios Scenario names (`"D0"`, `"D1"`, `"D2"`) and/or a named
#'   list of [contamination()] objects.
#' @param methods Named list of methods (see [default_methods()]).
#' @param times Evaluation grid (hours since waking).
#' @param references Optional named list of [sim_reference()] objects; the
#'   built-in reference is used for any trimester not supplied.
#'
#' @return An object of class `nmea_sim_design`.
#' @export
#' @examples
#' sim_design(trimesters = "T2", n = 50, seeds = 1:2, scenarios = c("D0", "D2"))
sim_design <- function(trimesters = c("T1", "T2", "T3"), n = 100, seeds = 1:10,
                       scenarios = c("D0", "D1", "D2"), methods = default_methods(),
                       times = seq(0, 18, by = 0.1), references = list()) {
  trimesters <- match.arg(trimesters, c("T1", "T2", "T3"), several.ok = TRUE)
  .check_scalar_number(n, "n", lower = 3)
  if (!is.numeric(times) || length(times) < 2 || any(!is.finite(times)) || any(diff(times) <= 0) ||
      times[1] < 0) {
    stop("`times` must be a finite, increasing, non-negative grid of at least two points.", call. = FALSE)
  }
  if (!is.list(seeds)) {
    if (!is.numeric(seeds) || any(!is.finite(seeds))) stop("`seeds` must be integers.", call. = FALSE)
    seeds <- lapply(stats::setNames(trimesters, trimesters), function(t) {
      as.integer((3 * as.double(seeds) + match(t, c("T1", "T2", "T3"))) %% 2147483647)
    })
  }
  if (!all(trimesters %in% names(seeds))) stop("`seeds` must have an entry for every trimester.", call. = FALSE)
  seeds <- lapply(seeds[trimesters], as.integer)
  if (anyDuplicated(unlist(seeds))) {
    warning("Some seeds are shared between trimesters; their cohorts are not independent.", call. = FALSE)
  }
  if (is.character(scenarios)) {
    scenarios <- stats::setNames(lapply(scenarios, .scenario_contamination), scenarios)
  }
  if (!is.list(scenarios) || is.null(names(scenarios)) ||
      !all(vapply(scenarios, inherits, logical(1), "nmea_contamination"))) {
    stop("`scenarios` must be scenario names or a named list of contamination() objects.", call. = FALSE)
  }
  if (!is.list(methods) || is.null(names(methods)) || anyDuplicated(names(methods)) ||
      !all(vapply(methods, inherits, logical(1), "nmea_sim_method"))) {
    stop("`methods` must be a uniquely named list of method specifications.", call. = FALSE)
  }
  refs <- lapply(stats::setNames(trimesters, trimesters), function(t) {
    if (!is.null(references[[t]])) references[[t]] else sim_reference(t)
  })
  structure(list(trimesters = trimesters, n = as.integer(n), seeds = seeds, scenarios = scenarios,
                 methods = methods, times = times, references = refs),
            class = "nmea_sim_design")
}

#' @export
print.nmea_sim_design <- function(x, ...) {
  n_cohorts <- sum(lengths(x$seeds))
  cat(sprintf("<nmea_sim_design> %s; %d subjects per cohort; %d base cohorts x %d scenarios = %d data sets\n",
              paste(x$trimesters, collapse = "/"), x$n, n_cohorts, length(x$scenarios),
              n_cohorts * length(x$scenarios)))
  cat("  scenarios:", paste(names(x$scenarios), collapse = ", "), "\n")
  cat("  methods  :", paste(names(x$methods), collapse = ", "), "\n")
  invisible(x)
}

#' Run a simulation study
#'
#' Generates every cohort and scenario of a [sim_design()], fits every method,
#' and keeps compact results (curves on the evaluation grid, measurement flags
#' and subject outcomes); full model objects are not kept. Failures are
#' recorded, never replaced.
#'
#' With `cache_dir`, each data set's results are saved as they finish and
#' reused when the same design is run again, so interrupted studies resume.
#' Cache entries are keyed by the design, the case and the package version.
#'
#' @param design An `nmea_sim_design`.
#' @param cache_dir Optional directory for per-case results.
#' @param parallel Logical; run base cohorts in parallel with
#'   `future.apply::future_lapply()` (configure workers with
#'   `future::plan()`). The package must be installed for workers to load it.
#' @param verbose Logical; print progress.
#'
#' @return An object of class `nmea_sim_runs` with the `design` and a list of
#'   `cases`, one per trimester, seed and scenario.
#' @export
sim_run <- function(design, cache_dir = NULL, parallel = FALSE, verbose = TRUE) {
  if (!inherits(design, "nmea_sim_design")) stop("`design` must come from sim_design().", call. = FALSE)
  jobs <- do.call(rbind, lapply(design$trimesters, function(t) {
    data.frame(trimester = t, seed = design$seeds[[t]], stringsAsFactors = FALSE)
  }))
  if (!is.null(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  key <- .design_key(design)
  run_job <- function(i) {
    .run_base_cohort(design, jobs$trimester[i], jobs$seed[i], cache_dir, key, verbose)
  }
  idx <- seq_len(nrow(jobs))
  cases <- if (isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package `future.apply` is required for parallel = TRUE.", call. = FALSE)
    }
    future.apply::future_lapply(idx, run_job, future.seed = NULL)
  } else {
    lapply(idx, run_job)
  }
  structure(list(design = design, cases = unlist(cases, recursive = FALSE)), class = "nmea_sim_runs")
}

.run_base_cohort <- function(design, trimester, seed, cache_dir, key, verbose) {
  base <- NULL
  out <- list()
  for (sc in names(design$scenarios)) {
    case_id <- sprintf("%s_%s_%s", trimester, seed, sc)
    path <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(case_id, "_", key, ".rds"))
    if (!is.null(path) && file.exists(path)) {
      out[[case_id]] <- readRDS(path)
      next
    }
    if (is.null(base)) base <- sim_cohort(design$references[[trimester]], n = design$n, seed = seed)
    sim <- sim_contaminate(base, design$scenarios[[sc]])
    truth_curves <- nmea_curve(design$times, sim$truth$psi)
    methods <- lapply(design$methods, .run_method, sim = sim, times = design$times)
    case <- list(id = case_id, trimester = trimester, seed = seed, scenario = sc,
                 truth_curves = truth_curves, points = sim$truth$points,
                 observed_time = sim$observed$time, methods = methods)
    if (!is.null(path)) {
      tmp <- paste0(path, ".tmp")
      saveRDS(case, tmp)
      file.rename(tmp, path)
    }
    if (isTRUE(verbose)) {
      status <- vapply(methods, `[[`, character(1), "status")
      message(sprintf("%s done (%s)", case_id, paste(names(status), status, sep = ": ", collapse = ", ")))
    }
    out[[case_id]] <- case
  }
  out
}

# Deparsed source of every function in the package namespace, so that cached
# results are not reused after the code changes (the development version
# number does not change with every edit).
.code_signature <- function() {
  ns <- asNamespace("cortisol.nmea")
  nms <- sort(ls(ns, all.names = TRUE))
  lapply(stats::setNames(nms, nms), function(n) {
    obj <- get(n, envir = ns)
    if (is.function(obj)) deparse(obj) else NULL
  })
}

.design_key <- function(design) {
  tmp <- tempfile()
  on.exit(unlink(tmp))
  # Methods hold closures; hash their plain-data description instead.
  method_sig <- lapply(design$methods, function(m) {
    if (m$type == "gamm") return(list("gamm", m$k))
    s <- m$steps
    list("nmea", if (is.null(s$outlier)) NULL else .spec_plain(s$outlier$spec), s$initial_alpha, s$min_obs,
         s$fvu_max, s$c1_positive, s$alpha_grid, unclass(m$control))
  })
  saveRDS(list(design$n, lapply(design$scenarios, unclass), method_sig, design$times,
               lapply(design$references, unclass),
               as.character(utils::packageVersion("cortisol.nmea")), .code_signature(),
               as.character(utils::packageVersion("saemix")),
               as.character(utils::packageVersion("mgcv"))), tmp)
  substr(unname(tools::md5sum(tmp)), 1, 10)
}

#' @export
print.nmea_sim_runs <- function(x, ...) {
  status <- unlist(lapply(x$cases, function(cs) vapply(cs$methods, `[[`, character(1), "status")))
  cat(sprintf("<nmea_sim_runs> %d data sets; %d method fits, %d failed\n",
              length(x$cases), length(status), sum(status != "ok")))
  print(x$design)
  invisible(x)
}
