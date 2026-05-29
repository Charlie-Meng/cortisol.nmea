test_that("GAM baseline fits and predicts on a time grid", {
  skip_if_not_installed("mgcv")

  data(toy_cortisol)
  toy_small <- subset(toy_cortisol, Subject %in% unique(Subject)[1:5])
  time_grid <- seq(0, 18, by = 3)
  subject_ids <- unique(as.character(toy_small$Subject))

  fit <- fit_gam_baseline(
    toy_small,
    id = "Subject",
    time = "Time",
    response = "Conc_Obs",
    k = 4
  )
  pred <- predict_gam_baseline(fit, subject_ids = subject_ids, time_grid = time_grid)

  expect_s3_class(fit, "nmea_gam_baseline")
  expect_equal(dim(pred), c(length(subject_ids), length(time_grid)))
  expect_true(all(is.finite(pred)))
})
