test_that("package citation is available", {
  cit <- utils::citation("cortisol.nmea")
  expect_s3_class(cit, "citation")
  expect_match(format(cit, style = "text"), "contributed equally")
})
