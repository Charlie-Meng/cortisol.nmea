# Synthetic-only package data generator.
# Run from the package root with:
#   Rscript data-raw/make_toy_cortisol.R

source("R/utils.R")
source("R/time_transform.R")
source("R/model_curve.R")
source("R/features.R")
source("R/simulate_parameters.R")
source("R/simulate_observed.R")

set.seed(20260529)

theta <- c(mu = 0.9, s = 1.1, c1 = 32, c0 = 2.2)
covariance <- matrix(
  c(
    0.18^2, 0.01,   0.12, 0.00,
    0.01,   0.15^2, 0.20, 0.01,
    0.12,   0.20,   5.50, 0.15,
    0.00,   0.01,   0.15, 0.45
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(c("mu", "s", "c1", "c0"), c("mu", "s", "c1", "c0"))
)

psi <- simulate_nmea_parameters(
  n = 30,
  theta = theta,
  covariance = covariance,
  alpha = 1.2,
  seed = 20260529
)

time_templates <- list(
  standard = c(0, 0.5, 1, 3, 6, 9, 12, 15),
  early_dense = c(0, 0.33, 0.75, 1.5, 4, 8, 12, 16),
  sparse = c(0, 0.75, 3, 8, 14)
)

sigma_fun <- function(time) 0.7 + 0.15 * sqrt(pmax(time, 0)) + 0.35 * exp(-((time - 1) / 0.8)^2)

toy_cortisol <- simulate_observed_data(
  psi = psi,
  time_templates = time_templates,
  sigma = sigma_fun,
  seed = 20260529
)

save(toy_cortisol, file = "data/toy_cortisol.rda", compress = "xz")
utils::write.csv(toy_cortisol, file = "inst/extdata/toy_cortisol.csv", row.names = FALSE)
