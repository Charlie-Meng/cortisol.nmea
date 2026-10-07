# Small synthetic cohort for tests. Parameters are arbitrary but plausible;
# no real participant data are involved.
make_test_cohort <- function(n = 24, seed = 1, spikes = 3) {
  set.seed(seed)
  psi <- cbind(
    mu = stats::rnorm(n, 0.5, 0.35),
    s = exp(stats::rnorm(n, log(0.75), 0.2)),
    c1 = exp(stats::rnorm(n, log(30), 0.3)),
    c0 = exp(stats::rnorm(n, log(2.5), 0.25)),
    alpha = 1.5
  )
  rownames(psi) <- sprintf("S%02d", seq_len(n))
  base_times <- c(0, 0.6, 3, 8, 13)
  rows <- lapply(seq_len(n), function(i) {
    tt <- sort(pmax(0, base_times + c(0, stats::runif(4, -0.2, 0.4))))
    truth <- nmea_curve(tt, psi[i, ])
    data.frame(subject = rownames(psi)[i], time = tt,
               y = pmax(0.01, truth + stats::rnorm(length(tt), 0, 1.2)), truth = truth)
  })
  d <- do.call(rbind, rows)
  spike_rows <- integer()
  if (spikes > 0) {
    spike_rows <- sample(nrow(d), spikes)
    d$y[spike_rows] <- d$y[spike_rows] + 15
  }
  list(data = cortisol_data(d, "subject", "time", "y"), psi = psi, spike_rows = spike_rows)
}

fast_control <- function() nmea_control(iterations = c(60, 30))
