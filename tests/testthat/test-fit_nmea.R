test_that("fit_nmea returns fitted parameters and observed-row predictions", {
  skip_if_not_installed("saemix")

  data(toy_cortisol)
  toy_small <- subset(toy_cortisol, Subject %in% unique(Subject)[1:4])

  fit <- fit_nmea(
    toy_small,
    id = "Subject",
    time = "Time",
    response = "Conc_Obs",
    alpha = 1.2,
    psi0 = c(mu = 1, s = 1, c1 = 30, c0 = 2),
    control = list(
      nbiter.saemix = c(20, 10),
      nbiter.burn = 10,
      nbiter.map = 10,
      fim = FALSE
    )
  )

  expect_s3_class(fit, "nmea_fit")
  expect_equal(nrow(fit$psi), 4)
  expect_true(all(c("mu", "s", "c1", "c0", "alpha") %in% colnames(fit$psi)))
  expect_length(fit$ypred, nrow(toy_small))
  expect_true(all(is.finite(fit$ypred)))
  expect_true(is.finite(fit$RSS))
})
