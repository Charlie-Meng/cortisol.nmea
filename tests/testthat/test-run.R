small_runs <- function(cache_dir = NULL) {
  ctrl <- fast_control()
  design <- sim_design(
    trimesters = c("T2", "T3"), n = 20, seeds = list(T2 = 1:2, T3 = 3:4), scenarios = c("D0", "D2"),
    methods = list(Direct = method_nmea(nmea_steps_direct(alpha_grid = c(1, 1.5)), ctrl),
                   Full = method_nmea(nmea_steps(alpha_grid = c(1, 1.5)), ctrl),
                   GAMM = method_gamm())
  )
  sim_run(design, cache_dir = cache_dir, verbose = FALSE)
}

test_that("sim_design validates its inputs", {
  d <- sim_design(trimesters = "T2", n = 10, seeds = 1:3, scenarios = c("D0", "D2"))
  expect_s3_class(d, "nmea_sim_design")
  expect_equal(d$seeds$T2, as.integer(3 * (1:3) + 2))
  expect_error(sim_design(scenarios = "D9"), "Unknown scenario")
  expect_error(sim_design(methods = list(method_gamm())), "named")
  custom <- sim_design(trimesters = "T1", seeds = 1, scenarios = list(neg = contamination(0.1, c(2, 3))))
  expect_identical(names(custom$scenarios), "neg")
})

test_that("a small study runs, caches and evaluates consistently", {
  skip_on_cran()
  cache <- withr::local_tempdir()
  runs <- small_runs(cache)
  expect_s3_class(runs, "nmea_sim_runs")
  expect_length(runs$cases, 8)
  expect_length(list.files(cache, pattern = "[.]rds$"), 8)

  again <- small_runs(cache)
  expect_identical(again$cases, runs$cases)

  subj <- eval_subjects(runs)
  expect_equal(nrow(subj), 8 * 3 * 20)
  expect_true(all(is.na(subj$rmse[!subj$available])))

  s <- eval_summary(runs)
  expect_true(all(c("rmse", "auc_error", "common_n") %in% s$summary$measure))
  expect_true(all(s$cohorts$common_n <= 20))
  pooled <- s$summary[s$summary$trimester == "All" & s$summary$measure == "rmse", ]
  expect_equal(nrow(pooled), 2 * 3)

  p <- eval_paired(runs, "Full", "Direct")
  expect_true(all(c("difference") %in% p$summary$measure))

  det <- eval_detection(runs, "Full")$cohorts
  expect_true(any(det$valid))
  d0 <- det[det$scenario == "D0" & det$valid, ]
  expect_true(all(d0$TP == 0 & d0$FN == 0))
  expect_error(eval_detection(runs, "GAMM"), "no outlier step")
  expect_error(eval_detection(runs, "Direct"), "no outlier step")

  ret <- eval_retention(runs, "Full")$cohorts
  pct <- ret[, c("retained", "excluded_min_obs", "excluded_fvu", "excluded_c1", "not_fitted")]
  expect_equal(unname(rowSums(pct)), rep(100, nrow(ret)))

  tn <- tune_threshold(runs, k = c(2, 3.5, 6), sd = "reference")
  d2 <- tn$summary[tn$summary$scenario == "D2", ]
  expect_true(all(diff(d2$FPR[d2$trimester == "T2"]) <= 0))

  expect_s3_class(plot_cohort_curves(runs, "T2", 1), "ggplot")
  expect_s3_class(plot_subject(runs, "T2_S0001", "T2", 1, "D2"), "ggplot")
  expect_s3_class(plot_recovery(runs), "ggplot")
  expect_s3_class(plot_detection(runs), "ggplot")
  expect_s3_class(plot_retention(runs), "ggplot")
})

test_that("the full-method flags in a run equal a direct pipeline call", {
  skip_on_cran()
  runs <- small_runs()
  cs <- runs$cases[["T3_4_D2"]]
  sim <- sim_contaminate(sim_cohort(sim_reference("T3"), n = 20, seed = 4), "D2")
  r <- nmea_pipeline(sim$observed, steps = nmea_steps(alpha_grid = c(1, 1.5)), control = fast_control())
  expect_identical(cs$methods$Full$points$flagged, r$points$flagged)
  expect_equal(cs$methods$Full$curves, predict(r, times = runs$design$times))
})
