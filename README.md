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
- Optional model fitting through `saemix`.
- Synthetic parameter, latent-curve, and observed-data simulation.
- Simulation V3-style missingness, sparsity, contamination, and cleaning
  utilities.
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
psi <- simulate_nmea_parameters(n = 20, seed = 11)

observed <- simulate_observed_data(
  psi,
  time_templates = list(
    standard = c(0, 0.5, 1, 3, 6, 9, 12, 15),
    sparse = c(0, 0.75, 3, 8, 14)
  ),
  sigma = function(time) 0.7 + 0.15 * sqrt(pmax(time, 0)),
  seed = 12
)

sparse <- apply_missing_scenario(
  observed,
  scenario = "S4_SparseMorningMissing",
  seed = 13
)

contaminated <- apply_contamination_scenario(
  sparse,
  scenario = "SC2_TimingError",
  seed = 14
)
```

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
