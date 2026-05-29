# cortisol.nmea

`cortisol.nmea` provides tools for nonlinear mixed-effects modeling of
diurnal cortisol curves and synthetic simulation benchmarks.

The package is designed as the public, reusable method layer for the cortisol
NMEA project. Private real-data analysis scripts and participant data should
remain outside this repository.

## What is included

- NMEA curve evaluation with interpretable parameters: `mu`, `s`, `c1`, `c0`,
  and `alpha`.
- Derived cortisol features: AUC, EML, PCL, AR, and DDC.
- Optional model fitting through `saemix`, including symmetric and alpha-grid
  wrappers.
- Optional smooth GAM baseline through `mgcv`.
- Synthetic parameter, latent-curve, and observed-data simulation.
- Simulation V3-style missingness, sparsity, contamination, cleaning, method
  comparison, and benchmark table utilities.
- A small fully synthetic toy data set for examples and vignettes.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Charlie-Meng/cortisol.nmea")
```

## Quick example

```r
library(cortisol.nmea)

psi <- simulate_nmea_parameters(n = 3, seed = 1)
time_grid <- seq(0, 18, by = 0.25)
curves <- model_eval(time_grid, psi)

features <- summary_features(psi)
features
```

## Synthetic simulation example

```r
study <- run_simulation_study(
  n = 20,
  missing_scenarios = c("S0_CleanObserved", "S4_SparseMorningMissing"),
  contamination_scenarios = c("C1_Spike", "SC2_TimingError"),
  seed = 11
)

comparison <- run_method_comparison(
  study,
  methods = c("OracleTruth", "SymmetricTruth")
)

comparison$summary
make_benchmark_tables(comparison$summary)
```

Lower-level functions such as `simulate_observed_data()`,
`apply_missing_scenario()`, and `apply_contamination_scenario()` remain
available for custom simulation designs.

## Package workflow status

The public API now covers curve evaluation, feature extraction, optional NMEA
fitting, synthetic Simulation V3-style data generation, method comparison, GAM
baseline fitting, and benchmark table construction. The current public package
intentionally excludes the private real-data analysis layer.

## Vignettes

- `vignette("nmea-model", package = "cortisol.nmea")`
- `vignette("simulation-workflow", package = "cortisol.nmea")`
- `vignette("benchmark-example", package = "cortisol.nmea")`

## Private data policy

This repository is public. Do not commit real raw data, private fitted
workspaces, real-data `.rda/.rds` files, derived private CSV/XLSX outputs,
PDFs, slide decks, or real-data analysis folders.

The public package should use only synthetic examples, toy data, and functions
that are safe to share.
