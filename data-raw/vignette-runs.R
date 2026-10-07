# Builds inst/extdata/vignette-runs.rds, the small precomputed study used by the
# simulation and threshold vignettes. Synthetic data only; takes a few minutes.
# Run from the package root: Rscript data-raw/vignette-runs.R
devtools::load_all(quiet = TRUE)
design <- sim_design(
  trimesters = "T2", n = 40, seeds = list(T2 = c(11L, 12L, 13L)),
  scenarios = c("D0", "D2"), methods = default_methods(),
  times = seq(0, 18, by = 0.25)
)
runs <- sim_run(design)
saveRDS(runs, "inst/extdata/vignette-runs.rds", compress = "xz")
cat("size:", file.info("inst/extdata/vignette-runs.rds")$size, "bytes\n")
