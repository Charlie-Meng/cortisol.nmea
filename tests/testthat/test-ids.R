# Regression tests from the PR #1 cross-review: subject identifiers must be
# matched as identifiers, whatever their form.

relabel <- function(co, f) {
  d <- as.data.frame(unclass(co$data), stringsAsFactors = FALSE)
  d$subject <- f(d$subject)
  d$point_id <- NULL
  cortisol_data(d, "subject", "time", "y")
}

test_that("numeric identifiers select subjects by identifier, not position", {
  skip_on_cran()
  co <- make_test_cohort(n = 20, spikes = 0)
  for (offset in c(0, 100)) {
    cd <- relabel(co, function(s) as.integer(sub("S", "", s)) + offset)
    fit <- nmea_fit(cd, control = fast_control())
    id <- 2 + offset
    by_number <- predict(fit, times = c(0, 1, 6), subjects = id)
    by_text <- predict(fit, times = c(0, 1, 6), subjects = as.character(id))
    expect_identical(rownames(by_number), as.character(id))
    expect_identical(by_number, by_text)
    expect_equal(unname(by_number[1, ]), nmea_curve(c(0, 1, 6), coef(fit)[as.character(id), ]))
    expect_error(predict(fit, times = 0, subjects = 999), "Unknown subject")
  }
})

test_that("nmea_result prediction matches identifiers and rejects unknown ones", {
  skip_on_cran()
  co <- make_test_cohort(n = 20, spikes = 0)
  cd <- relabel(co, function(s) as.integer(sub("S", "", s)))
  res <- nmea_pipeline(cd, steps = nmea_steps_direct(alpha_grid = 1), control = fast_control())
  expect_identical(predict(res, times = 0, subjects = 2), predict(res, times = 0, subjects = "2"))
  expect_error(predict(res, times = 0, subjects = "nobody"), "Unknown subject")
})

test_that("GAMM handles identifiers containing separators without changing the fit", {
  skip_on_cran()
  co <- make_test_cohort(n = 20, spikes = 0)
  plain <- fit_gamm_sanchez(co$data, times = c(0, 2, 8))
  odd <- fit_gamm_sanchez(relabel(co, function(s) paste0("site/", s)), times = c(0, 2, 8))
  expect_identical(odd$status, "ok")
  expect_identical(rownames(odd$curves), paste0("site/", rownames(plain$curves)))
  expect_equal(unname(odd$curves), unname(plain$curves))
  # Different prefixes, same suffix: must stay distinct subjects.
  mixed <- relabel(co, function(s) ifelse(as.integer(sub("S", "", s)) %% 2 == 0,
                                          paste0("a/", s), paste0("b/", s)))
  expect_identical(fit_gamm_sanchez(mixed, times = 0)$status, "ok")
})

test_that("alpha search keeps diagnostics for every candidate", {
  skip_on_cran()
  co <- make_test_cohort(n = 20, spikes = 0)
  a <- nmea_fit_alpha(co$data, alpha_grid = c(1, 1.5), control = fast_control())
  expect_true(all(c("fim_ok", "n_warnings", "warnings", "error", "elapsed") %in% names(a$profile)))
  expect_equal(nrow(a$profile), 2)
  expect_null(a$fits)
})
