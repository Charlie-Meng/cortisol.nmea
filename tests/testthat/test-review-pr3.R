# Regression tests from the PR #3 cross-review.

one_design <- function(methods, references = list(), scenarios = "D2", n = 20) {
  sim_design(trimesters = "T2", n = n, seeds = list(T2 = 7L), scenarios = scenarios,
             methods = methods, references = references)
}

test_that("a custom rule without key is never served from the cache (fix review 2026-10-28)", {
  skip_on_cran()
  cfg <- new.env()
  cfg$cutoff <- 3
  rule <- outlier_custom(function(fit) abs(fit$residual) > cfg$cutoff)
  d <- one_design(list(A = method_nmea(nmea_steps(outlier = rule, alpha_grid = 1), fast_control())))
  cache <- withr::local_tempdir()
  expect_warning(r3 <- sim_run(d, cache_dir = cache, verbose = FALSE), "Caching is off")
  expect_length(list.files(cache), 0)
  cfg$cutoff <- 9
  expect_warning(r9 <- sim_run(d, cache_dir = cache, verbose = FALSE), "Caching is off")
  fresh <- suppressWarnings(sim_run(d, verbose = FALSE))
  expect_identical(r9$cases[[1]]$methods$A$points$flagged, fresh$cases[[1]]$methods$A$points$flagged)
  expect_lt(sum(r9$cases[[1]]$methods$A$points$flagged), sum(r3$cases[[1]]$methods$A$points$flagged))
})

test_that("a failed iterative refit keeps its diagnostics (fix review 2026-10-28)", {
  skip_on_cran()
  sim <- sim_contaminate(sim_cohort(sim_reference("T2"), n = 20, seed = 7), "D2")
  d <- sim$observed
  d$trimester <- NULL                    # unlabelled data -> iterative SD
  real_fit <- nmea_fit
  calls <- 0
  local_mocked_bindings(nmea_fit = function(...) {
    calls <<- calls + 1
    if (calls == 2) {                    # the iterative SD refit
      return(structure(list(status = "failure", error = "injected", warnings = "w1",
                            elapsed = 1.5), class = "nmea_fit"))
    }
    real_fit(...)
  })
  res <- nmea_pipeline(d, steps = nmea_steps(alpha_grid = 1), control = fast_control())
  expect_identical(res$status, "failure")
  expect_match(res$error, "iterative SD refit failed: injected")
  expect_true(all(is.na(res$points$flagged)))
  expect_identical(res$stages$rule$status, "failure")
  expect_identical(res$stages$rule$warnings, "w1")
  expect_equal(res$stages$rule$elapsed, 1.5)
  expect_true(res$stages$rule$n_kept > 0)
})

test_that("a design's own reference is used both for fitting and for tuning", {
  skip_on_cran()
  noisy <- sim_reference_modify(sim_reference("T2"), noise_multiplier = 2)
  ctrl <- fast_control()
  d <- one_design(list(Full = method_nmea(nmea_steps(alpha_grid = 1), ctrl)), list(T2 = noisy))
  runs <- sim_run(d, verbose = FALSE)
  cs <- runs$cases[[1]]
  expect_equal(cs$methods$Full$cutoff, 3.5 * noisy$noise_sd)
  tn <- tune_threshold(runs, k = 3.5, sd = "reference")$cohorts
  expect_equal(tn$cutoff, cs$methods$Full$cutoff)
  expect_equal(tn$flagged, sum(cs$methods$Full$points$flagged))
})

test_that("a failing outlier rule becomes a failure of that method only", {
  skip_on_cran()
  ctrl <- fast_control()
  boom <- outlier_sd(3, function(fit) stop("injected calibration failure"))
  d <- one_design(list(Full = method_nmea(nmea_steps(outlier = boom, alpha_grid = 1), ctrl),
                       Direct = method_nmea(nmea_steps_direct(alpha_grid = 1), ctrl)))
  runs <- sim_run(d, verbose = FALSE)
  full <- runs$cases[[1]]$methods$Full
  expect_identical(full$status, "failure")
  expect_match(full$error, "injected calibration failure")
  expect_true(all(is.na(full$points$flagged)))
  expect_identical(runs$cases[[1]]$methods$Direct$status, "ok")
  # Invalid for detection, shown as "not evaluated", and plots still work.
  expect_false(eval_detection(runs, "Full")$cohorts$valid)
  p <- plot_subject(runs, "T2_S0001", "T2", 7, "D2")
  expect_true(all(p$layers[[2]]$data$class == "NE"))
  expect_s3_class(plot_cohort_curves(runs, "T2", 7, methods = "Full"), "ggplot")
})

test_that("iterative SD diagnostics are kept", {
  skip_on_cran()
  ctrl <- fast_control()
  d <- one_design(list(Full = method_nmea(nmea_steps(outlier = outlier_sd(3, "iterative"),
                                                     alpha_grid = 1), ctrl)))
  diag <- sim_run(d, verbose = FALSE)$cases[[1]]$methods$Full$diagnostics$rule
  expect_identical(diag$stage, "iterative_sd")
  expect_true(all(c("status", "fim_ok", "warnings", "elapsed", "n_kept") %in% names(diag)))
})

test_that("failed initial fits leave tuning and cohort plots usable", {
  skip_on_cran()
  bad <- nmea_control(iterations = c(5, 5), psi0 = c(mu = 1, s = 0, c1 = 1, c0 = 0))
  d <- one_design(list(Full = method_nmea(nmea_steps(alpha_grid = 1), bad)))
  runs <- sim_run(d, verbose = FALSE)
  expect_identical(runs$cases[[1]]$methods$Full$status, "failure")
  tn <- tune_threshold(runs, k = c(3, 3.5))
  expect_false(any(tn$cohorts$valid))
  expect_equal(tn$summary$valid, c(0L, 0L))
  expect_equal(tn$summary$planned, c(1L, 1L))
  expect_s3_class(plot_cohort_curves(runs, "T2", 7), "ggplot")
})

test_that("invalid evaluation grids are rejected", {
  expect_error(sim_design(times = c(0, 2, 1)), "increasing")
  expect_error(sim_design(times = c(0, NA)), "increasing")
})
