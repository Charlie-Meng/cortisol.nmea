test_that("references are valid and contain aggregates only", {
  for (tri in c("T1", "T2", "T3")) {
    ref <- sim_reference(tri)
    expect_s3_class(ref, "nmea_reference")
    expect_equal(dim(ref$knots), c(4L, 11L))
    expect_equal(sum(ref$patterns), 1, tolerance = 1e-4)
    expect_gte(ref$n_reference, 100)
  }
  expect_equal(sim_reference("T2")$alpha, 1.5)
  expect_equal(sim_reference("T1")$alpha, 1.3)
})

test_that("sim_reference_modify changes components and validates", {
  ref <- sim_reference_modify(sim_reference("T2"), alpha = 1.2, noise_multiplier = 2)
  expect_equal(ref$alpha, 1.2)
  expect_equal(ref$noise_sd, 2 * sim_reference("T2")$noise_sd)
  expect_error(sim_reference_modify(sim_reference("T2"), nonsense = 1), "Unknown")
  expect_error(sim_reference_modify(sim_reference("T2"), alpha = -1), "Invalid")
})

test_that("sim_cohort is reproducible and keeps truth separate from observations", {
  a <- sim_cohort(sim_reference("T2"), n = 30, seed = 11)
  b <- sim_cohort(sim_reference("T2"), n = 30, seed = 11)
  expect_identical(a$observed, b$observed)
  expect_s3_class(a$observed, "cortisol_data")
  expect_false(any(c("truth", "contaminated") %in% names(a$observed)))
  expect_equal(nrow(a$truth$psi), 30)
  expect_true(all(a$truth$psi[, "alpha"] == 1.5))
  expect_true(all(a$observed$y >= 0.01))
  expect_true(all(tapply(a$observed$time, a$observed$subject, min) == 0))
  expect_true(all(tapply(a$observed$time, a$observed$subject, function(t) all(diff(t) > 0))))
  expect_equal(a$truth$points$truth, unname(unlist(lapply(rownames(a$truth$psi), function(id) {
    nmea_curve(a$observed$time[a$observed$subject == id], a$truth$psi[id, ])
  }))))
})

test_that("simulation does not change the user's random number stream", {
  set.seed(5)
  expected <- stats::runif(1)
  set.seed(5)
  sim_contaminate(sim_cohort(sim_reference("T1"), n = 5, seed = 1), "D2")
  expect_identical(stats::runif(1), expected)
})

test_that("D0, D1 and D2 share the base cohort and follow their definitions", {
  base <- sim_cohort(sim_reference("T3"), n = 100, seed = 3)
  d0 <- sim_contaminate(base, "D0")
  d1 <- sim_contaminate(base, "D1")
  d2 <- sim_contaminate(base, "D2")
  expect_identical(d0$observed$y, base$observed$y)
  for (d in list(d1, d2)) {
    expect_identical(d$observed$time, base$observed$time)
    expect_identical(d$truth$points$clean_y, base$truth$points$clean_y)
  }
  p2 <- d2$truth$points
  expect_equal(sum(p2$contaminated), round(0.10 * nrow(p2)))
  ratio <- p2$delta[p2$contaminated] / p2$sigma[p2$contaminated]
  expect_true(all(ratio >= 8 - 1e-9 & ratio <= 12 + 1e-9))
  expect_true(all(d2$observed$y[!p2$contaminated] == base$observed$y[!p2$contaminated]))
  r1 <- d1$truth$points$delta[d1$truth$points$contaminated] / d1$truth$points$sigma[d1$truth$points$contaminated]
  expect_true(all(r1 >= 4 - 1e-9 & r1 <= 8 + 1e-9))
})

test_that("custom and timing contamination work", {
  base <- sim_cohort(sim_reference("T2"), n = 40, seed = 4)
  neg <- sim_contaminate(base, contamination(0.1, c(1, 2), scale = "absolute", direction = "negative"))
  expect_true(all(neg$truth$points$delta <= 0))
  win <- sim_contaminate(base, contamination(0.05, c(5, 6), time_window = c(0, 1)))
  expect_true(all(win$truth$points$time[win$truth$points$contaminated] <= 1))
  tm <- sim_contaminate(base, contamination(0.1, c(0.5, 1.5), type = "timing"))
  moved <- tm$truth$points$contaminated
  expect_true(any(moved))
  expect_true(all(tm$observed$time >= 0))
  expect_identical(tm$observed$y, base$observed$y)
  expect_false(any(tm$observed$time[moved] == base$observed$time[moved]))
})

test_that("sim_thin keeps at least min_obs per subject", {
  base <- sim_cohort(sim_reference("T2"), n = 30, seed = 8)
  th <- sim_thin(base, mcar = 0.5, morning = 0.3, min_obs = 3)
  expect_true(all(table(th$observed$subject) >= 3))
  expect_identical(th$truth$points$point_id, th$observed$point_id)
  sp <- sim_thin(base, keep = 3)
  expect_true(all(table(sp$observed$subject) == 3))
})
