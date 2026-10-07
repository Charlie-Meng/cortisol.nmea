# cortisol.nmea

<!-- badges: start -->
[![R-CMD-check](https://github.com/Charlie-Meng/cortisol.nmea/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Charlie-Meng/cortisol.nmea/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Salivary cortisol is usually collected only four or five times a day, at
irregular times. `cortisol.nmea` estimates each person's full diurnal curve from
such sparse data with an asymmetric scaled-logistic nonlinear mixed-effects
model (NMEA), and provides the simulation framework we used to evaluate it.

The package has two uses:

- **Fit your own data.** Run the NMEA workflow on long-format cortisol data:
  an initial fit, removal of aberrant measurements, subject screening and a
  search over the asymmetry of the morning rise. Every step can be switched off
  and every tuning parameter changed. The fitted curves give the cortisol
  awakening level, peak, awakening rise, decline slope and area under the curve.
- **Run simulation studies.** Generate synthetic pregnancy cohorts with known
  true curves, add contamination, compare methods and summarize curve recovery,
  outlier detection and subject retention, with figures and an HTML report.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Charlie-Meng/cortisol.nmea")
```

## Fitting your own data

```r
library(cortisol.nmea)

d <- cortisol_data(my_data, id = "participant", time = "hours_since_waking",
                   value = "cortisol_nmol_l")
fit <- nmea_pipeline(d)
fit                          # what was removed or excluded, and why
curves <- predict(fit)       # curves on a 0-18 h grid
nmea_features(fit)           # AUC, EML, PCL, AR, DDC per subject
```

Time is hours since waking (the first sample of the day), concentrations are in
nmol/L (`convert_cortisol_units()` converts from ug/dL or ng/mL). Fit each
trimester or study wave separately. See `vignette("cortisol-nmea")`.

## Simulation studies

```r
design <- sim_design(trimesters = c("T1", "T2", "T3"), n = 100, seeds = 1:10,
                     scenarios = c("D0", "D1", "D2"))
runs <- sim_run(design, cache_dir = "sim-cache")
eval_summary(runs)$summary
plot_recovery(runs)
sim_report(runs, "simulation-report.html")
```

See `vignette("simulation")`, `vignette("thresholds")` and
`vignette("paper-simulation")`.

## Data

The package contains no participant data. Simulated cohorts are generated from
aggregate summaries of a prenatal cortisol cohort: quantiles and rank
correlations of the individual curve parameters, the sampling schedule and
missingness, and the measurement noise level. A privacy check in continuous
integration blocks data files from entering the repository.

## Comparison method

`fit_gamm_sanchez()` adapts the generalized additive mixed model of Sánchez et
al. (2012), *Modeling the salivary cortisol profile in population research: the
Multi-Ethnic Study of Atherosclerosis*, American Journal of Epidemiology
176(10):918-928, <https://doi.org/10.1093/aje/kws182>. It is our implementation
for comparison, not the authors' code.

## Citation

Use `citation("cortisol.nmea")`. Yuntian Meng and Yu Gu contributed equally;
Xing Qiu supervised the work.

## License

MIT © Yuntian Meng, Yu Gu and Xing Qiu
