test_that("benchmark table helpers produce winners and wide metric tables", {
  summary_df <- data.frame(
    Scenario = rep(c("S0", "S4"), each = 2),
    ScenarioGroup = rep(c("Clean", "Sparse"), each = 2),
    Method = rep(c("A", "B"), times = 2),
    Curve_RMSE = c(0.1, 0.2, 0.5, 0.4),
    AUC_MAE = c(1, 0.8, 2.5, 3),
    stringsAsFactors = FALSE
  )

  winners <- summarize_winners(summary_df)
  expect_equal(winners$Winner_Curve_RMSE, c("A", "B"))
  expect_equal(winners$Winner_AUC_MAE, c("B", "A"))

  tables <- make_benchmark_tables(summary_df)
  expect_named(tables, c("long", "wide", "winners"))
  expect_equal(nrow(tables$long), 8)
  expect_true(any(grepl("^Value\\.", names(tables$wide))))
})
