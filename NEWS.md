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
  uncleaned fit. Outlier rules: `outlier_fixed()`, `outlier_sd()`,
  `outlier_custom()`. The default rule is provisionally `outlier_fixed(6)`
  until the `k x sd` default is chosen.
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
