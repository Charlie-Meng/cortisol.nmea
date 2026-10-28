psi <- c(mu = 0.6, s = 0.75, c1 = 30, c0 = 2.5, alpha = 1.5)

test_that("time warping is invertible and leaves post-peak times unchanged", {
  tt <- seq(0, 3, by = 0.25)
  w <- nmea_warp_time(tt, mu = 0.6, alpha = 1.5)
  expect_equal(nmea_unwarp_time(w, mu = 0.6, alpha = 1.5), tt)
  expect_identical(w[tt >= 0.6], tt[tt >= 0.6])
})

test_that("stable curve equals the original formula where the latter is finite", {
  tt <- seq(0, 18, by = 0.1)
  m <- matrix(psi, nrow = 1, dimnames = list("a", names(psi)))
  original <- .nmea_saemix_model(m, rep(1L, length(tt)), matrix(tt, ncol = 1))
  expect_equal(nmea_curve(tt, psi), unname(original), tolerance = 1e-12)
})

test_that("stable curve does not overflow for extreme parameters", {
  extreme <- c(mu = 15, s = 0.001, c1 = 30, c0 = 2.5, alpha = 1)
  m <- matrix(extreme, nrow = 1, dimnames = list("a", names(extreme)))
  original <- suppressWarnings(.nmea_saemix_model(m, 1L, matrix(0, ncol = 1)))
  expect_true(is.nan(original))
  expect_equal(nmea_curve(0, extreme), 2.5)
})

test_that("nmea_curve returns one row per subject", {
  p <- rbind(a = psi, b = psi * c(1, 1, 2, 1, 1))
  out <- nmea_curve(c(0, 1, 2), p)
  expect_equal(dim(out), c(2L, 3L))
  expect_identical(rownames(out), c("a", "b"))
})
