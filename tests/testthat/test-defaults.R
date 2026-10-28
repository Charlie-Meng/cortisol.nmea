test_that("default rule uses the reference SD for trimester data and the iterative SD otherwise", {
  skip_on_cran()
  sim <- sim_contaminate(sim_cohort(sim_reference("T2"), n = 25, seed = 3), "D2")
  fit <- nmea_fit(sim$observed, control = fast_control())
  r <- .apply_outlier_rule(outlier_default(), fit)
  expect_equal(r$cutoff, 3.5 * sim_reference("T2")$noise_sd)

  unlabelled <- sim$observed
  unlabelled$trimester <- NULL
  fit2 <- nmea_fit(unlabelled, control = fast_control())
  keep <- abs(fit2$residual) <= 2.75 * stats::mad(fit2$residual)
  refit <- nmea_fit(.subset_cd(fit2$data, keep), alpha = 1, control = fast_control())
  expect_equal(.apply_outlier_rule(outlier_default(), fit2)$cutoff, 3 * refit$sigma)
  expect_error(.apply_outlier_rule(outlier_sd(3.5, "reference"), fit2), "trimester")
})

test_that("nmea_steps() defaults to outlier_default()", {
  expect_identical(nmea_steps()$outlier$spec$type, "default")
})

test_that("rules built on user functions are cacheable only with a key", {
  f <- function(fit) abs(fit$residual) > 3
  expect_false(.rule_cacheable(outlier_custom(f)))
  expect_true(.rule_cacheable(outlier_custom(f, key = "cut3")))
  expect_false(.rule_cacheable(outlier_sd(3, function(fit) 1)))
  expect_true(.rule_cacheable(outlier_sd(3, function(fit) 1, key = "sd1")))
  for (r in list(outlier_default(), outlier_fixed(6), outlier_sd(3.5, "reference"),
                 outlier_sd(3, "iterative"), outlier_sd(2.75, "mad"), NULL)) {
    expect_true(.rule_cacheable(r))
  }
  expect_error(outlier_custom(f, key = ""), "key")
  design_with <- function(key) {
    sim_design(trimesters = "T2", seeds = 1, scenarios = "D0",
               methods = list(A = method_nmea(nmea_steps(outlier = outlier_custom(f, key = key)))))
  }
  expect_false(identical(.design_key(design_with("a")), .design_key(design_with("b"))))
  expect_identical(.design_key(design_with("a")), .design_key(design_with("a")))
})

test_that("a seed vector gives independent trimester-specific seeds", {
  d <- sim_design(seeds = 1:10)
  s <- unlist(d$seeds)
  expect_false(anyDuplicated(s) > 0)
  expect_equal(d$seeds$T2, as.integer(3 * (1:10) + 2))
  expect_warning(sim_design(seeds = list(T1 = 1:2, T2 = 2:3, T3 = 5:6)), "not independent")
})

test_that("pooled summaries are NA when a planned trimester has no valid cohort", {
  expect_true(is.na(.pooled(c(1, 3, 5, 7, NA, NA), rep(c("T1", "T2", "T3"), each = 2))[["mean"]]))
  full <- .pooled(c(1, 3, 5, 7, 2, 4), rep(c("T1", "T2", "T3"), each = 2))
  expect_equal(full[["mean"]], mean(c(2, 6, 3)))
})
