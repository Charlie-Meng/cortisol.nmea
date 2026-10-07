# cortisol.nmea (development version)

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
