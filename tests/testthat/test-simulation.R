test_that("synthetic parameter and observed-data generators work", {
  psi <- simulate_nmea_parameters(n = 5, seed = 1)
  expect_equal(nrow(psi), 5)
  expect_true(all(psi[, "s"] > 0))
  expect_true(all(psi[, "c1"] > 0))

  obs <- simulate_observed_data(
    psi,
    time_templates = list(c(0, 0.5, 1, 3, 6)),
    sigma = 0.5,
    seed = 2
  )
  expect_equal(length(unique(obs$Subject)), 5)
  expect_true(all(c("Truth", "Sigma_t", "Conc_Obs") %in% names(obs)))
})

test_that("missingness and contamination scenarios add audit columns", {
  psi <- simulate_nmea_parameters(n = 4, seed = 3)
  obs <- simulate_observed_data(
    psi,
    time_templates = list(c(0, 0.5, 1, 3, 6, 9)),
    sigma = 0.5,
    seed = 4
  )

  sparse <- apply_missing_scenario(obs, scenario = "S4_SparseMorningMissing", seed = 5)
  expect_true(all(table(sparse$Subject) >= 3))
  expect_equal(unique(sparse$Scenario), "S4_SparseMorningMissing")

  contam <- apply_contamination_scenario(sparse, scenario = "SC2_TimingError", seed = 6)
  expect_true(all(c("Time_Clean", "Conc_Obs_Clean", "Contam_Class", "Contam_Flag") %in% names(contam)))
  expect_true(any(contam$Contam_Flag == 1))
})
