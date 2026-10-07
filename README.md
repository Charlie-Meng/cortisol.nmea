# cortisol.nmea

<!-- badges: start -->
[![R-CMD-check](https://github.com/Charlie-Meng/cortisol.nmea/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Charlie-Meng/cortisol.nmea/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> **Status:** under active rewrite on the `rewrite` branch. The API is not
> stable yet; the first release will be 0.1.0.

`cortisol.nmea` fits an asymmetric scaled-logistic nonlinear mixed-effects
model (NMEA) to sparse, irregularly sampled diurnal cortisol data and provides
a simulation framework for evaluating curve recovery under contamination.

## Planned functionality

- **Method layer:** bring your own cortisol data and run the NMEA workflow
  (initial fit, residual-based outlier removal with a `k x sd` threshold,
  minimum-observation and FVU screening, refit, asymmetry search). Every step
  can be switched off and every tuning parameter can be changed.
- **Simulation layer:** generate synthetic pregnancy-trimester cohorts from
  aggregate summaries, add configurable contamination, compare methods,
  evaluate recovery, detection and retention, and produce figures and reports.

## Comparison methods

The generalized additive mixed model comparator is an adaptation of the
approach of Sánchez et al. (2012), *American Journal of Epidemiology*
176(10):918-928, <doi:10.1093/aje/kws182>. It is not the original authors'
code.

## Data policy

This repository is public. It contains **no participant-level data**. The
simulation defaults are aggregate summaries only. Do not commit raw data,
fitted participant-level results, workspaces, or derived private outputs.
`tools/check_privacy.R` enforces this in continuous integration.

## Authors

Yuntian Meng (maintainer), Yu Gu, and Xing Qiu. Yuntian Meng and Yu Gu
contributed equally.

## License

MIT
