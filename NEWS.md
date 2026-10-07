# cortisol.nmea 0.0.0.9000

* Complete rewrite of the package around the final NMEA simulation workflow.
  The earlier 0.0.0.9000 prototype (May 2026) is preserved in the `master`
  branch history and is not compatible with this version.

## Method layer

* `cortisol_data()`, `convert_cortisol_units()` and `prep_cortisol()` validate
  and prepare long-format cortisol data.
* `nmea_curve()`, `nmea_warp_time()` and `nmea_unwarp_time()` evaluate the
  asymmetric scaled-logistic curve in a numerically stable form.
* `nmea_features()` returns AUC, EML, PCL, AR and DDC. The AUC is the exact
  piecewise integral of the warped curve; it corrects the original closed
  form, which was exact only when `0 < mu < tmax`.
* `nmea_fit()` and `nmea_fit_alpha()` fit the NMEA model with saemix; the
  saemix structural model is arithmetically identical to the original code so
  that published results are reproducible.
* `nmea_pipeline()` runs the configurable workflow (initial fit, outlier
  removal, minimum observations, FVU and `c1 > 0` screening, asymmetry search);
  `nmea_steps()` switches steps on or off and `nmea_steps_direct()` gives the
  uncleaned fit. Outlier rules: `outlier_default()`, `outlier_fixed()`,
  `outlier_sd()` (`sd` = number, `"reference"`, `"iterative"`, `"model"`,
  `"mad"` or a function) and `outlier_custom()`. The default rule uses
  `3.5 x` the trimester reference noise SD for T1/T2/T3 data and otherwise
  `3 x` an iterative SD estimate; `outlier_fixed(6)` reproduces the original
  analysis.
* `fit_gamm_sanchez()` provides the GAMM comparator adapted from Sánchez et al.
  (2012).
* Subject identifiers are always matched as identifiers: `predict()` no longer
  treats numeric identifiers as row positions, and the GAMM comparator accepts
  identifiers containing `/` (cross-review of PR #1).
* `nmea_fit_alpha()` keeps compact diagnostics (FIM availability, warnings,
  error, time) for every alpha candidate.

## Simulation layer

* `sim_reference()` returns aggregate pregnancy-trimester summaries (parameter
  quantile knots with a Gaussian copula, sampling protocol and missingness
  patterns, noise SD by time). No participant-level data are included.
  `sim_reference_modify()` customizes any component.
* `sim_cohort()` generates synthetic cohorts with truth kept separate from the
  observed data; `sim_contaminate()` adds the paper's D0/D1/D2 scenarios or a
  custom `contamination()` (spikes in either direction, absolute or
  noise-scaled, optional time window, or timing errors); `sim_thin()` adds
  optional missingness and sparsity.
* Cross-review of PR #2: contamination never raises a measurement through the
  floor (relevant to log-normal noise and custom floors); `corr` must be a
  symmetric correlation matrix with unit diagonal; the split of rare
  missingness patterns is a documented uniform assumption instead of an
  estimate from fewer than 20 people.

## Running, evaluating and reporting studies

* `sim_design()` and `sim_run()` run complete studies (trimesters x seeds x
  scenarios x methods) with optional on-disk caching for resuming and optional
  parallel execution; failures are recorded, never replaced.
  `method_nmea()`, `method_gamm()`, `default_methods()` and
  `methods_threshold_grid()` specify the compared methods.
* `eval_subjects()`, `eval_summary()` (common-subject curve RMSE and curve-AUC
  error with Monte Carlo intervals), `eval_paired()`, `eval_detection()` and
  `eval_retention()` summarize results.
* `tune_threshold()` reports detection trade-offs over `k` without refitting.
* `plot_cohort_curves()`, `plot_subject()`, `plot_recovery()`,
  `plot_detection()` and `plot_retention()` reproduce the paper's figure
  panels; `sim_report()` writes an HTML report.
* Cross-review of PR #3: cache keys now describe custom rules when the key is
  computed, including the global variables they read (codetools); rules using
  `sd = "reference"` take the simulation design's own reference, so fitting and
  `tune_threshold()` agree; a failing outlier rule (e.g. the iterative refit)
  becomes a recorded failure instead of stopping `sim_run()`, and its
  calibration diagnostics are kept; measurements without a removal decision
  are plotted as "not evaluated"; tuning, cohort plots and the report handle
  failed or single-method studies; summaries flag incomplete plans
  (`complete`); evaluation grids are validated.
