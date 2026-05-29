test_that("time transformation leaves post-peak times unchanged", {
  expect_equal(time_forward(c(0, 1, 2), alpha = 2, mu = 1), c(-1, 1, 2))
  expect_equal(time_backward(c(-1, 1, 2), alpha = 2, mu = 1), c(0, 1, 2))
})

test_that("NMEA curve evaluation returns one row per subject", {
  psi <- matrix(
    c(1, 1, 30, 2, 1.2, 1.3, 0.8, 25, 1.5, 1.1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("A", "B"), c("mu", "s", "c1", "c0", "alpha"))
  )
  pred <- model_eval(c(0, 1, 2), psi)
  expect_equal(dim(pred), c(2, 3))
  expect_true(all(is.finite(pred)))
  expect_equal(rownames(pred), c("A", "B"))
})

test_that("summary features are returned for each subject", {
  psi <- matrix(
    c(1, 1, 30, 2, 1.2),
    nrow = 1,
    dimnames = list("A", c("mu", "s", "c1", "c0", "alpha"))
  )
  features <- summary_features(psi)
  expect_equal(nrow(features), 1)
  expect_true(all(c("AUC", "EML", "PCL", "AR", "DDC") %in% names(features)))
})
