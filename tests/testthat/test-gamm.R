test_that("GAMM comparator predicts curves with a verified basis", {
  skip_on_cran()
  co <- make_test_cohort(n = 20, spikes = 0)
  g <- fit_gamm_sanchez(co$data, times = c(0, 1, 6, 12))
  expect_identical(g$status, "ok")
  expect_equal(dim(g$curves), c(20L, 4L))
  expect_lt(g$basis_check, 1e-6)
  expect_equal(length(g$fitted), nrow(co$data))
})
