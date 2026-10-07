test_that("cortisol_data standardizes columns and keeps row order", {
  raw <- data.frame(pid = c("b", "b", "b", "a", "a", "a"), h = c(0, 1, 5, 0, 2, 9),
                    v = c(10, 15, 4, 12, 9, 3))
  d <- cortisol_data(raw, "pid", "h", "v")
  expect_s3_class(d, "cortisol_data")
  expect_identical(d$subject, raw$pid)
  expect_identical(d$point_id, c("b_001", "b_002", "b_003", "a_001", "a_002", "a_003"))
})

test_that("cortisol_data rejects invalid input", {
  raw <- data.frame(id = rep(1:2, each = 3), t = c(0, 1, 2, 0, 1, 2), y = c(1, 2, 3, 4, 5, 6))
  expect_error(cortisol_data(raw, "id", "time", "y"), "not found")
  bad <- raw
  bad$t[2] <- -1
  expect_error(cortisol_data(bad, "id", "t", "y"), "non-negative")
  bad <- raw
  bad$y[1] <- NA
  expect_error(cortisol_data(bad, "id", "t", "y"), "finite")
  expect_warning(cortisol_data(raw[-1, ], "id", "t", "y"), "fewer than 3")
})

test_that("unit conversion uses the molar mass of cortisol", {
  expect_equal(convert_cortisol_units(1, "ug/dL"), 27.589, tolerance = 1e-4)
  expect_equal(convert_cortisol_units(1, "ng/mL"), 2.7589, tolerance = 1e-4)
  expect_identical(convert_cortisol_units(5, "nmol/L"), 5)
})

test_that("prep_cortisol applies the time window", {
  raw <- data.frame(id = rep(1:3, each = 3), t = rep(c(0, 8, 20), 3), y = 1:9)
  d <- suppressWarnings(prep_cortisol(cortisol_data(raw, "id", "t", "y")))
  expect_true(all(d$time <= 18))
  expect_equal(nrow(d), 6)
})
