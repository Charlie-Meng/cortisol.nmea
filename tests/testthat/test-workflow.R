test_that("high-level simulation and method comparison workflow returns tables", {
  study <- run_simulation_study(
    n = 8,
    missing_scenarios = c("S0_CleanObserved", "S4_SparseMorningMissing"),
    contamination_scenarios = c("C1_Spike", "SC2_TimingError"),
    seed = 7
  )

  expect_s3_class(study, "nmea_simulation")
  expect_true(all(c("S0_CleanObserved", "S4_SparseMorningMissing", "C1_Spike", "SC2_TimingError") %in% names(study$scenarios)))

  comparison <- run_method_comparison(
    study,
    methods = c("OracleTruth", "SymmetricTruth")
  )

  expect_s3_class(comparison, "nmea_method_comparison")
  expect_gt(nrow(comparison$summary), 0)
  expect_true(all(c("Scenario", "Method", "Curve_RMSE", "AUC_MAE", "Status") %in% names(comparison$summary)))

  winners <- summarize_winners(comparison$summary)
  expect_true("Winner_Curve_RMSE" %in% names(winners))

  tables <- make_benchmark_tables(comparison$summary)
  expect_named(tables, c("long", "wide", "winners"))
})
