fake_fit <- function(res) list(residual = res, sigma = 2, data = data.frame(y = res))

test_that("fixed and sd rules flag by absolute residual", {
  res <- c(-7, -1, 0, 5, 6.5)
  expect_identical(.apply_outlier_rule(outlier_fixed(6), fake_fit(res))$flag,
                   c(TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(.apply_outlier_rule(outlier_sd(k = 3, sd = "model"), fake_fit(res))$cutoff, 6)
  expect_equal(.apply_outlier_rule(outlier_sd(k = 3.5, sd = 1.5), fake_fit(res))$cutoff, 5.25)
  expect_equal(.apply_outlier_rule(outlier_sd(2, "mad"), fake_fit(res))$cutoff,
               2 * stats::mad(res))
})

test_that("numbers are treated as fixed cutoffs and invalid rules are rejected", {
  expect_identical(.apply_outlier_rule(6, fake_fit(c(7, 1)))$flag, c(TRUE, FALSE))
  expect_error(outlier_sd(3.5), "must be supplied")
  expect_error(outlier_sd(3.5, "iqr"), "must be a number")
  bad <- outlier_custom(function(fit) c(TRUE, NA))
  expect_error(.apply_outlier_rule(bad, fake_fit(c(1, 2))), "logical")
})
