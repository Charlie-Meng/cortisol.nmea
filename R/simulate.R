#' Simulation reference for a pregnancy trimester
#'
#' Returns the aggregate summaries used to generate synthetic cohorts. The
#' defaults were computed from a prenatal cortisol cohort and contain **no
#' participant-level values**:
#'
#' * `knots`, `knot_probs`, `corr`: individual parameters (`mu`, `log s`,
#'   `log c1`, `log c0`) as 11 marginal quantile knots (2.5% to 97.5%) joined by
#'   a Gaussian copula with normal-score correlation `corr`. Marginal tails
#'   beyond the outer knots are extrapolated linearly.
#' * `alpha`: the asymmetry of the generated (true) curves.
#' * `protocol`, `deviation_knots`, `patterns`: the sampling design. Sample `j`
#'   is taken at `protocol[j]` hours after waking plus a deviation drawn from
#'   the quantile knots in row `j` of `deviation_knots` (clamped to the
#'   2.5%-97.5% range); `patterns` gives the probability of each set of
#'   collected samples (e.g. `"1234"`: the bedtime sample is missing).
#' * `noise_time`, `noise_sigma`, `noise_sd`: the residual noise SD
#'   \eqn{\sigma(t)}, linearly interpolated between knots and constant beyond
#'   them, and its overall level.
#' * `floor`: lower limit applied to simulated concentrations.
#'
#' The reference parameters are post-quality-control empirical Bayes
#' estimates, so they describe screened subjects with shrinkage, not an
#' unscreened physiological population.
#'
#' Any component can be changed with [sim_reference_modify()].
#'
#' @param trimester `"T1"`, `"T2"` or `"T3"`.
#'
#' @return An object of class `nmea_reference`.
#' @export
#' @examples
#' ref <- sim_reference("T2")
#' ref
sim_reference <- function(trimester = c("T1", "T2", "T3")) {
  trimester <- match.arg(trimester)
  d <- .nmea_reference_defaults[[trimester]]
  rownames(d$corr) <- colnames(d$corr) <- rownames(d$knots)
  d$floor <- 0.01
  d$source <- "aggregate summaries of a prenatal cortisol cohort (post-QC EBEs)"
  .validate_reference(structure(d, class = "nmea_reference"))
}

#' Modify a simulation reference
#'
#' Changes components of an [sim_reference()] object, e.g. to simulate a
#' population with a different asymmetry, a noisier assay, or a different
#' sampling protocol.
#'
#' @param reference An `nmea_reference`.
#' @param ... Named components to replace (see [sim_reference()] for names).
#'   `noise_multiplier = m` is a shortcut that scales `noise_sigma` and
#'   `noise_sd` by `m`.
#'
#' @return The modified `nmea_reference`.
#' @export
#' @examples
#' ref <- sim_reference_modify(sim_reference("T2"), alpha = 1.2, noise_multiplier = 1.5)
#' ref$noise_sd
sim_reference_modify <- function(reference, ...) {
  if (!inherits(reference, "nmea_reference")) stop("`reference` must be an nmea_reference.", call. = FALSE)
  changes <- list(...)
  if (!is.null(changes$noise_multiplier)) {
    m <- changes$noise_multiplier
    .check_scalar_number(m, "noise_multiplier", lower = 0)
    reference$noise_sigma <- reference$noise_sigma * m
    reference$noise_sd <- reference$noise_sd * m
    changes$noise_multiplier <- NULL
  }
  unknown <- setdiff(names(changes), names(reference))
  if (length(unknown)) stop("Unknown reference component(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  for (nm in names(changes)) reference[[nm]] <- changes[[nm]]
  reference$source <- paste(reference$source, "(modified)")
  .validate_reference(reference)
}

.validate_reference <- function(ref) {
  ok <- is.numeric(ref$alpha) && length(ref$alpha) == 1 && ref$alpha > 0 &&
    is.matrix(ref$knots) && nrow(ref$knots) == 4 && ncol(ref$knots) == length(ref$knot_probs) &&
    all(diff(ref$knot_probs) > 0) && all(apply(ref$knots, 1, function(k) all(diff(k) >= 0))) &&
    is.matrix(ref$corr) && all(dim(ref$corr) == 4) &&
    is.numeric(ref$protocol) && nrow(ref$deviation_knots) == length(ref$protocol) &&
    ncol(ref$deviation_knots) == length(ref$knot_probs) &&
    !is.null(names(ref$patterns)) && abs(sum(ref$patterns) - 1) < 1e-4 && all(ref$patterns >= 0) &&
    length(ref$noise_time) == length(ref$noise_sigma) && all(ref$noise_sigma > 0) &&
    is.numeric(ref$floor) && length(ref$floor) == 1
  if (!ok) stop("Invalid simulation reference.", call. = FALSE)
  idx <- unique(as.integer(unlist(strsplit(names(ref$patterns), ""))))
  if (any(is.na(idx)) || any(idx < 1 | idx > length(ref$protocol))) {
    stop("Sampling patterns must list sample numbers present in `protocol`.", call. = FALSE)
  }
  if (inherits(try(chol(ref$corr), silent = TRUE), "try-error")) {
    stop("`corr` must be a positive definite correlation matrix.", call. = FALSE)
  }
  ref
}

#' @export
print.nmea_reference <- function(x, ...) {
  cat(sprintf("<nmea_reference> %s: alpha = %s, noise SD = %.3f nmol/L\n",
              x$name, format(x$alpha), x$noise_sd))
  med <- x$knots[, which.min(abs(x$knot_probs - 0.5))]
  cat(sprintf("  median parameters: mu = %.2f h, s = %.2f, c1 = %.1f, c0 = %.2f\n",
              med[1], exp(med[2]), exp(med[3]), exp(med[4])))
  cat("  protocol (h):", paste(x$protocol, collapse = ", "), "\n")
  cat(sprintf("  all samples collected: %.0f%%\n", 100 * sum(x$patterns[nchar(names(x$patterns)) == length(x$protocol)])))
  cat("  source:", x$source, "\n")
  invisible(x)
}

# Inverse CDF from quantile knots; linear interpolation, optional linear tails.
.knot_quantile <- function(u, probs, knots, extrapolate = TRUE) {
  k <- length(probs)
  x <- stats::approx(probs, knots, xout = pmin(pmax(u, probs[1]), probs[k]), ties = "ordered")$y
  if (extrapolate) {
    lo <- u < probs[1]
    hi <- u > probs[k]
    x[lo] <- knots[1] - (knots[2] - knots[1]) / (probs[2] - probs[1]) * (probs[1] - u[lo])
    x[hi] <- knots[k] + (knots[k] - knots[k - 1]) / (probs[k] - probs[k - 1]) * (u[hi] - probs[k])
  }
  x
}

.noise_sd_at <- function(reference, t) {
  stats::approx(reference$noise_time, reference$noise_sigma, xout = t, rule = 2, ties = mean)$y
}

#' Simulate a synthetic cortisol cohort
#'
#' Draws individual NMEA parameters, sampling times and measurement noise from
#' a simulation reference. The observed data and the truth are stored
#' separately; fitting functions should only ever receive `$observed`.
#'
#' Random numbers are drawn in a fixed order (parameters, sampling design,
#' noise) from `seed`; the user's global random number stream is left
#' unchanged.
#'
#' @param reference An `nmea_reference` from [sim_reference()].
#' @param n Number of subjects.
#' @param seed Integer seed.
#' @param noise `"additive"` (default): \eqn{y = \max(floor, f(t) + \sigma(t) Z)},
#'   which slightly biases low concentrations upward through the floor;
#'   `"lognormal"`: a positive, mean-preserving multiplicative alternative with
#'   the same variance.
#'
#' @return An object of class `cortisol_sim` with elements `observed` (a
#'   [cortisol_data()] object), `truth` (list with `psi`, the true parameter
#'   matrix, and `points`, one row per measurement with the true curve value,
#'   noise SD, clean value and contamination fields), `reference`, `seed`,
#'   `scenario` and `contamination`.
#' @export
#' @examples
#' sim <- sim_cohort(sim_reference("T2"), n = 5, seed = 1)
#' sim$observed
sim_cohort <- function(reference, n, seed, noise = c("additive", "lognormal")) {
  if (!inherits(reference, "nmea_reference")) stop("`reference` must come from sim_reference().", call. = FALSE)
  .check_scalar_number(n, "n", lower = 1)
  .check_scalar_number(seed, "seed")
  noise <- match.arg(noise)
  n <- as.integer(n)
  .with_preserved_rng({
    set.seed(seed)
    psi <- .draw_parameters(reference, n)
    design <- .draw_design(reference, n)
    .assemble_cohort(reference, psi, design, noise, seed)
  })
}

.draw_parameters <- function(reference, n) {
  z <- matrix(stats::rnorm(n * 4), n, 4) %*% chol(reference$corr)
  u <- stats::pnorm(z)
  x <- vapply(1:4, function(j) .knot_quantile(u[, j], reference$knot_probs, reference$knots[j, ]),
              numeric(n))
  x <- matrix(x, nrow = n)
  psi <- cbind(mu = x[, 1], s = exp(x[, 2]), c1 = exp(x[, 3]), c0 = exp(x[, 4]),
               alpha = reference$alpha)
  rownames(psi) <- sprintf("%s_S%04d", reference$name, seq_len(n))
  psi
}

.draw_design <- function(reference, n) {
  pattern <- sample(names(reference$patterns), n, replace = TRUE, prob = reference$patterns)
  lapply(pattern, function(p) {
    idx <- sort(as.integer(strsplit(p, "")[[1]]))
    for (attempt in 1:100) {
      dev <- vapply(idx, function(j) {
        .knot_quantile(stats::runif(1), reference$knot_probs, reference$deviation_knots[j, ],
                       extrapolate = FALSE)
      }, numeric(1))
      tt <- pmax(0, reference$protocol[idx] + dev)
      tt <- tt - min(tt)   # time since the first collected sample (taken as waking)
      if (all(diff(tt) > 0)) return(list(sample = idx, time = round(tt, 4)))
    }
    stop("Could not draw increasing sampling times.", call. = FALSE)
  })
}

.assemble_cohort <- function(reference, psi, design, noise, seed) {
  ids <- rownames(psi)
  rows <- lapply(seq_along(ids), function(i) {
    tt <- design[[i]]$time
    truth <- nmea_curve(tt, psi[i, ])
    data.frame(subject = ids[i], sample = design[[i]]$sample, time = tt, truth = truth,
               stringsAsFactors = FALSE)
  })
  pts <- do.call(rbind, rows)
  pts$point_id <- paste0(pts$subject, "_", pts$sample)
  pts$sigma <- .noise_sd_at(reference, pts$time)
  innovation <- stats::rnorm(nrow(pts))
  if (noise == "additive") {
    raw <- pts$truth + pts$sigma * innovation
    pts$clean_y <- pmax(reference$floor, raw)
    pts$floor_applied <- raw < reference$floor
  } else {
    lv <- log1p((pts$sigma / pts$truth)^2)
    pts$clean_y <- pts$truth * exp(sqrt(lv) * innovation - lv / 2)
    pts$floor_applied <- FALSE
  }
  pts$y <- pts$clean_y
  pts$contaminated <- FALSE
  pts$contamination_type <- "none"
  pts$delta <- 0
  pts$trimester <- reference$name
  observed <- cortisol_data(pts, "subject", "time", "y", point_id = "point_id",
                            trimester = "trimester", sample = "sample")
  structure(list(
    observed = observed,
    truth = list(psi = psi, points = pts[, c("point_id", "subject", "sample", "time", "truth",
                                              "sigma", "clean_y", "floor_applied", "y",
                                              "contaminated", "contamination_type", "delta")]),
    reference = reference, seed = seed, noise = noise, scenario = "D0", contamination = NULL
  ), class = "cortisol_sim")
}

#' @export
print.cortisol_sim <- function(x, ...) {
  pts <- x$truth$points
  cat(sprintf("<cortisol_sim> %s, scenario %s: %d subjects, %d measurements, seed %s\n",
              x$reference$name, x$scenario, nrow(x$truth$psi), nrow(pts), format(x$seed)))
  cat(sprintf("  contaminated measurements: %d (%.1f%%); floor applied: %.1f%%\n",
              sum(pts$contaminated), 100 * mean(pts$contaminated), 100 * mean(pts$floor_applied)))
  invisible(x)
}

#' Contamination settings
#'
#' Describes measurement contamination added by [sim_contaminate()].
#'
#' * `type = "spike"` adds `delta` to a random `rate` fraction of all
#'   measurements. `delta` is drawn uniformly from `magnitude`; with
#'   `scale = "sigma"` it is multiplied by the noise SD \eqn{\sigma(t)} at that
#'   time, with `scale = "absolute"` it is in nmol/L. `direction` controls the
#'   sign. This is the contamination used in the paper.
#' * `type = "timing"` shifts the recorded time of the selected measurements by
#'   a uniform `magnitude` hours (random sign, never before waking), mimicking
#'   misreported sampling times. Removing such points does not repair them.
#'
#' `time_window` restricts which measurements are eligible.
#'
#' @param rate Fraction of measurements to contaminate (at least one).
#' @param magnitude Length-two range of the contamination size.
#' @param type `"spike"` or `"timing"`.
#' @param scale `"sigma"` or `"absolute"` (spikes only).
#' @param direction `"positive"`, `"negative"` or `"both"` (spikes only).
#' @param time_window Optional length-two range of eligible times (hours).
#'
#' @return An object of class `nmea_contamination`.
#' @export
#' @examples
#' contamination(rate = 0.10, magnitude = c(8, 12))   # the paper's D2
#' contamination(rate = 0.05, magnitude = c(0.5, 1.5), type = "timing")
contamination <- function(rate, magnitude, type = c("spike", "timing"),
                          scale = c("sigma", "absolute"),
                          direction = c("positive", "negative", "both"),
                          time_window = NULL) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  direction <- match.arg(direction)
  .check_scalar_number(rate, "rate", lower = 0, upper = 1)
  if (!is.numeric(magnitude) || length(magnitude) != 2L || any(!is.finite(magnitude)) ||
      magnitude[1] > magnitude[2] || magnitude[1] < 0) {
    stop("`magnitude` must be a non-negative increasing range of length 2.", call. = FALSE)
  }
  if (!is.null(time_window) && (!is.numeric(time_window) || length(time_window) != 2L)) {
    stop("`time_window` must be NULL or a numeric range of length 2.", call. = FALSE)
  }
  structure(list(rate = rate, magnitude = magnitude, type = type, scale = scale,
                 direction = direction, time_window = time_window),
            class = "nmea_contamination")
}

#' @export
print.nmea_contamination <- function(x, ...) {
  if (x$rate == 0) {
    cat("<nmea_contamination> none\n")
  } else if (x$type == "spike") {
    cat(sprintf("<nmea_contamination> %s spikes on %.1f%% of measurements, size U(%s, %s) x %s\n",
                x$direction, 100 * x$rate, x$magnitude[1], x$magnitude[2],
                if (x$scale == "sigma") "sigma(t)" else "nmol/L"))
  } else {
    cat(sprintf("<nmea_contamination> timing shifts on %.1f%% of measurements, U(%s, %s) h\n",
                100 * x$rate, x$magnitude[1], x$magnitude[2]))
  }
  invisible(x)
}

.scenario_contamination <- function(scenario) {
  switch(scenario,
    D0 = contamination(rate = 0, magnitude = c(0, 0)),
    D1 = contamination(rate = 0.05, magnitude = c(4, 8)),
    D2 = contamination(rate = 0.10, magnitude = c(8, 12)),
    stop("Unknown scenario; use \"D0\", \"D1\", \"D2\" or contamination().", call. = FALSE)
  )
}

#' Add contamination to a simulated cohort
#'
#' Contaminates the clean measurements of a [sim_cohort()] result. Scenarios
#' derived from the same cohort share subjects, sampling times and clean noise,
#' so methods can be compared pairwise across scenarios.
#'
#' Predefined scenarios: `"D0"` none; `"D1"` 5% positive spikes of
#' U(4, 8) x \eqn{\sigma(t)}; `"D2"` 10% positive spikes of U(8, 12) x
#' \eqn{\sigma(t)}.
#'
#' @param sim A `cortisol_sim` object (contamination always starts from its
#'   clean values).
#' @param scenario `"D0"`, `"D1"`, `"D2"`, or an object from [contamination()].
#' @param seed Integer seed; by default `sim$seed + 100`.
#'
#' @return A `cortisol_sim` with contaminated `observed` data and per-measurement
#'   contamination labels in `truth$points`.
#' @export
#' @examples
#' base <- sim_cohort(sim_reference("T2"), n = 10, seed = 1)
#' d2 <- sim_contaminate(base, "D2")
#' d2
sim_contaminate <- function(sim, scenario = "D2", seed = sim$seed + 100) {
  if (!inherits(sim, "cortisol_sim")) stop("`sim` must come from sim_cohort().", call. = FALSE)
  spec <- if (inherits(scenario, "nmea_contamination")) scenario else .scenario_contamination(scenario)
  .check_scalar_number(seed, "seed")
  ref <- sim$reference
  pts <- sim$truth$points
  pts$y <- pts$clean_y
  pts$contaminated <- FALSE
  pts$contamination_type <- "none"
  pts$delta <- 0
  time <- sim$truth$points$time
  # Recorded times start from the clean design.
  base_times <- pts$time

  if (spec$rate > 0) {
    .with_preserved_rng({
      set.seed(seed)
      eligible <- seq_len(nrow(pts))
      if (!is.null(spec$time_window)) {
        eligible <- eligible[pts$time >= spec$time_window[1] & pts$time <= spec$time_window[2]]
      }
      if (!length(eligible)) stop("No measurements fall in `time_window`.", call. = FALSE)
      k <- min(length(eligible), max(1L, round(spec$rate * nrow(pts))))
      j <- eligible[sample.int(length(eligible), k)]
      size <- stats::runif(k, spec$magnitude[1], spec$magnitude[2])
      if (spec$type == "spike") {
        sign <- switch(spec$direction, positive = rep(1, k), negative = rep(-1, k),
                       both = sample(c(-1, 1), k, replace = TRUE))
        delta <- sign * size * if (spec$scale == "sigma") pts$sigma[j] else 1
        pts$y[j] <- pmax(ref$floor, pts$clean_y[j] + delta)
        pts$delta[j] <- pts$y[j] - pts$clean_y[j]
        pts$contamination_type[j] <- "spike"
      } else {
        sign <- sample(c(-1, 1), k, replace = TRUE)
        shifted <- base_times[j] + sign * size
        flip <- shifted < 0
        shifted[flip] <- base_times[j][flip] + size[flip]
        time[j] <- shifted
        pts$delta[j] <- shifted - base_times[j]
        pts$contamination_type[j] <- "timing"
      }
    })
    pts$contaminated <- abs(pts$delta) > 1e-10
    pts$contamination_type[!pts$contaminated] <- "none"
  }

  obs <- data.frame(subject = pts$subject, time = time, y = pts$y, point_id = pts$point_id,
                    trimester = ref$name, sample = pts$sample, stringsAsFactors = FALSE)
  sim$observed <- cortisol_data(obs, "subject", "time", "y", point_id = "point_id",
                                trimester = "trimester", sample = "sample")
  sim$truth$points <- pts
  sim$scenario <- if (is.character(scenario)) scenario else "custom"
  sim$contamination <- spec
  sim
}

#' Remove measurements from a simulated cohort
#'
#' Optional sparsity and missingness settings from earlier simulation
#' experiments (not used in the paper's main design). Removal is applied to
#' both the observed data and the truth table and keeps at least `min_obs`
#' measurements per subject.
#'
#' @param sim A `cortisol_sim` object.
#' @param mcar Probability that any measurement is missing completely at random.
#' @param morning Additional missingness probability for measurements within
#'   the first hour after waking.
#' @param keep Optional number of measurements to keep per subject (chosen at
#'   random), e.g. `3` for very sparse sampling.
#' @param min_obs Minimum number of measurements kept per subject.
#' @param seed Integer seed; by default `sim$seed + 200`.
#'
#' @return A thinned `cortisol_sim`.
#' @export
#' @examples
#' base <- sim_cohort(sim_reference("T2"), n = 10, seed = 1)
#' sim_thin(base, mcar = 0.15)
sim_thin <- function(sim, mcar = 0, morning = 0, keep = NULL, min_obs = 3L, seed = sim$seed + 200) {
  if (!inherits(sim, "cortisol_sim")) stop("`sim` must come from sim_cohort().", call. = FALSE)
  .check_scalar_number(mcar, "mcar", 0, 1)
  .check_scalar_number(morning, "morning", 0, 1)
  .check_scalar_number(min_obs, "min_obs", lower = 1)
  pts <- sim$truth$points
  obs <- sim$observed
  kept <- .with_preserved_rng({
    set.seed(seed)
    unlist(lapply(split(seq_len(nrow(obs)), obs$subject), function(i) {
      p <- mcar + ifelse(obs$time[i] < 1, morning, 0)
      drop <- stats::runif(length(i)) < pmin(p, 1)
      cand <- i[!drop]
      if (!is.null(keep) && length(cand) > keep) cand <- sort(cand[sample.int(length(cand), keep)])
      if (length(cand) < min_obs) {
        extra <- setdiff(i, cand)
        need <- min(min_obs, length(i)) - length(cand)
        cand <- sort(c(cand, extra[sample.int(length(extra), need)]))
      }
      cand
    }), use.names = FALSE)
  })
  kept <- sort(kept)
  sim$observed <- .subset_cd(obs, kept)
  sim$truth$points <- pts[match(sim$observed$point_id, pts$point_id), , drop = FALSE]
  rownames(sim$truth$points) <- NULL
  sim
}
