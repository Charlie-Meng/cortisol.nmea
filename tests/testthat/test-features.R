numeric_auc <- function(p, tmax = 18) {
  stats::integrate(function(t) nmea_curve(t, p), 0, tmax, rel.tol = 1e-10,
                   subdivisions = 1000)$value
}

test_that("exact AUC matches numerical integration for any peak location", {
  for (mu in c(-0.4, 0.6, 19)) {
    p <- c(mu = mu, s = 0.8, c1 = 28, c0 = 2.2, alpha = 1.4)
    expect_equal(unname(nmea_features(p)[, "AUC"]), numeric_auc(p), tolerance = 1e-7)
  }
})

test_that("EML, PCL, AR and DDC agree with the curve", {
  p <- c(mu = 0.6, s = 0.8, c1 = 28, c0 = 2.2, alpha = 1.4)
  f <- nmea_features(p)
  expect_equal(unname(f[, "EML"]), nmea_curve(0, p))
  expect_equal(unname(f[, "PCL"]), nmea_curve(0.6, p))
  slope <- diff(nmea_curve(seq(0.6, 4, by = 1e-4), p)) / 1e-4
  expect_equal(unname(f[, "DDC"]), min(slope), tolerance = 1e-4)
  q <- p
  q["mu"] <- -0.3
  expect_equal(unname(nmea_features(q)[, "AR"]), 0)
})
