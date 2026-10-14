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

test_that("closures with different captured values get different cache keys", {
  make_rule <- function(cutoff) outlier_custom(function(fit) abs(fit$residual) > cutoff)
  d3 <- sim_design(trimesters = "T2", seeds = 1, scenarios = "D0",
                   methods = list(A = method_nmea(nmea_steps(outlier = make_rule(3)))))
  d9 <- sim_design(trimesters = "T2", seeds = 1, scenarios = "D0",
                   methods = list(A = method_nmea(nmea_steps(outlier = make_rule(9)))))
  expect_false(identical(.design_key(d3), .design_key(d9)))
  sd_rule <- function(v) outlier_sd(3, function(fit) v)
  e1 <- sim_design(trimesters = "T2", seeds = 1, scenarios = "D0",
                   methods = list(A = method_nmea(nmea_steps(outlier = sd_rule(1)))))
  e2 <- sim_design(trimesters = "T2", seeds = 1, scenarios = "D0",
                   methods = list(A = method_nmea(nmea_steps(outlier = sd_rule(2)))))
  expect_false(identical(.design_key(e1), .design_key(e2)))
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
