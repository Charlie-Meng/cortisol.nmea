test_that("nmea_fit returns ordered data and finite parameters", {
  skip_on_cran()
  co <- make_test_cohort(n = 20, spikes = 0)
  fit <- nmea_fit(co$data, alpha = 1, control = fast_control())
  expect_identical(fit$status, "ok")
  expect_identical(rownames(fit$psi), sort(unique(co$data$subject)))
  expect_true(all(is.finite(fit$psi)))
  expect_equal(length(fit$fvu), 20)
})

test_that("fitting does not change the user's random number stream", {
  skip_on_cran()
  co <- make_test_cohort(n = 20, spikes = 0)
  set.seed(99)
  expected <- stats::runif(1)
  set.seed(99)
  nmea_fit(co$data, control = fast_control())
  observed <- stats::runif(1)
  expect_identical(observed, expected)
})

test_that("full workflow records flags and subject outcomes", {
  skip_on_cran()
  co <- make_test_cohort(n = 24, spikes = 3)
  res <- nmea_pipeline(co$data, steps = nmea_steps(alpha_grid = c(1, 1.5)),
                       control = fast_control())
  expect_identical(res$status, "ok")
  expect_equal(nrow(res$points), nrow(co$data))
  expect_true(all(res$points$flagged == (abs(res$points$initial_residual) > res$cutoff)))
  expect_setequal(res$subjects$subject[res$subjects$outcome == "retained"], rownames(coef(res)))
  expect_true(all(res$subjects$outcome %in%
                    c("retained", "excluded_min_obs", "excluded_fvu", "excluded_c1")))
  curves <- predict(res, times = c(0, 1, 2))
  expect_equal(dim(curves), c(24L, 3L))
  expect_true(all(is.na(curves[res$subjects$outcome != "retained", ])))
})

test_that("direct workflow keeps every measurement and subject", {
  skip_on_cran()
  co <- make_test_cohort(n = 20, spikes = 2)
  res <- nmea_pipeline(co$data, steps = nmea_steps_direct(alpha_grid = c(1, 1.5)),
                       control = fast_control())
  expect_identical(res$status, "ok")
  expect_false(any(res$points$flagged))
  expect_true(all(res$subjects$outcome == "retained"))
  expect_null(res$stages$initial)
  expect_null(res$stages$screening)
})

test_that("too few subjects is reported as a failure, not an error", {
  skip_on_cran()
  co <- make_test_cohort(n = 4, spikes = 0)
  res <- nmea_pipeline(co$data, steps = nmea_steps(outlier = 0.01, alpha_grid = 1),
                       control = fast_control())
  expect_identical(res$status, "failure")
  expect_match(res$error, "Fewer than 3 subjects")
})
