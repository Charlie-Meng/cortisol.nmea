#' Apply missingness or sparse-sampling scenarios
#'
#' Public, generic version of the Simulation V3 missing/sparse generators.
#'
#' @param data Synthetic observed data with columns `Subject`, `Time`, and
#'   `Conc_Obs`.
#' @param scenario One of `S0_CleanObserved`, `S1_MCAR15`,
#'   `S2_MorningMissing`, `S3_Sparse34`, or `S4_SparseMorningMissing`.
#' @param seed Optional random seed.
#' @param min_keep Minimum observations retained per subject.
#' @param early_window Early morning window used by morning-missing scenarios.
#'
#' @return A data frame with retained observations and scenario labels.
#' @export
apply_missing_scenario <- function(data,
                                   scenario = c(
                                     "S0_CleanObserved",
                                     "S1_MCAR15",
                                     "S2_MorningMissing",
                                     "S3_Sparse34",
                                     "S4_SparseMorningMissing"
                                   ),
                                   seed = NULL,
                                   min_keep = 3,
                                   early_window = c(0, 1.5)) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)
  data <- as.data.frame(data)
  if (!all(c("Subject", "Time", "Conc_Obs") %in% names(data))) {
    stop("`data` must contain Subject, Time, and Conc_Obs columns.", call. = FALSE)
  }

  apply_one <- function(dat_sub) {
    n <- nrow(dat_sub)
    keep <- rep(TRUE, n)

    if (scenario %in% c("S1_MCAR15")) {
      keep <- stats::runif(n) > 0.15
    } else if (scenario %in% c("S2_MorningMissing")) {
      p_drop <- ifelse(
        dat_sub$Time >= early_window[1] & dat_sub$Time <= early_window[2],
        0.55,
        0.10
      )
      keep <- stats::runif(n) > p_drop
    } else if (scenario %in% c("S3_Sparse34")) {
      target_n <- min(n, sample(c(3L, 4L), size = 1))
      keep <- seq_len(n) %in% sort(sample(seq_len(n), size = target_n))
    } else if (scenario %in% c("S4_SparseMorningMissing")) {
      p_drop <- ifelse(
        dat_sub$Time >= early_window[1] & dat_sub$Time <= early_window[2],
        0.55,
        0.10
      )
      keep <- stats::runif(n) > p_drop
      keep <- .ensure_min_keep(keep, min_keep = min_keep)
      idx <- which(keep)
      target_n <- min(length(idx), sample(c(3L, 4L), size = 1))
      keep2 <- rep(FALSE, n)
      keep2[sort(sample(idx, size = target_n))] <- TRUE
      keep <- keep2
    }

    keep <- .ensure_min_keep(keep, min_keep = min_keep)
    dat_sub[keep, , drop = FALSE]
  }

  out <- do.call(rbind, lapply(split(data, data$Subject), apply_one))
  rownames(out) <- NULL
  out$Scenario <- scenario
  out$ScenarioGroup <- "MissingSparse"
  out
}

#' Apply synthetic contamination scenarios
#'
#' Public, generic version of the Simulation V3 contamination generator.
#'
#' @param data Synthetic observed data with `Subject`, `Time`, `Conc_Obs`, and
#'   preferably `Sigma_t`.
#' @param scenario One of `C0_None`, `C1_Spike`, `C2_TimingError`, `C3_Mixed`,
#'   `SC0_None`, `SC1_Spike`, `SC2_TimingError`, or `SC3_Mixed`.
#' @param seed Optional random seed.
#' @param spike_frac Fraction of rows affected by measurement spikes.
#' @param timing_frac Fraction of rows affected by timing errors.
#' @param time_range Allowed time range after timing shifts.
#' @param early_window Early window where timing errors are sampled more often.
#'
#' @return A data frame with contamination audit columns.
#' @export
apply_contamination_scenario <- function(data,
                                         scenario = c(
                                           "C0_None", "C1_Spike", "C2_TimingError", "C3_Mixed",
                                           "SC0_None", "SC1_Spike", "SC2_TimingError", "SC3_Mixed"
                                         ),
                                         seed = NULL,
                                         spike_frac = 0.05,
                                         timing_frac = 0.08,
                                         time_range = c(0, 18),
                                         early_window = c(0, 1.5)) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)
  data <- as.data.frame(data)
  if (!all(c("Subject", "Time", "Conc_Obs") %in% names(data))) {
    stop("`data` must contain Subject, Time, and Conc_Obs columns.", call. = FALSE)
  }
  if (!"Sigma_t" %in% names(data)) data$Sigma_t <- stats::sd(data$Conc_Obs, na.rm = TRUE)
  if (!is.finite(data$Sigma_t[1])) data$Sigma_t <- 1

  out <- data
  out$Time_Clean <- out$Time
  out$Conc_Obs_Clean <- out$Conc_Obs
  out$Contam_Class <- "none"
  out$Contam_Flag <- 0L

  add_spike <- function(dat) {
    n <- nrow(dat)
    m <- max(1L, floor(spike_frac * n))
    idx <- sample(seq_len(n), size = m)
    mult <- stats::runif(m, min = 2.5, max = 4.5)
    direction <- ifelse(stats::runif(m) <= 0.65, 1, -1)
    bump <- direction * mult * dat$Sigma_t[idx]
    dat$Conc_Obs[idx] <- pmax(dat$Conc_Obs[idx] + bump, 0.01)
    dat$Contam_Class[idx] <- "measurement_spike"
    dat$Contam_Flag[idx] <- 1L
    dat
  }

  add_timing <- function(dat) {
    n <- nrow(dat)
    m <- max(1L, floor(timing_frac * n))
    weight <- ifelse(dat$Time >= early_window[1] & dat$Time <= early_window[2], 3, 1)
    idx <- sample(seq_len(n), size = m, prob = weight)
    shift_abs <- ifelse(
      dat$Time[idx] >= early_window[1] & dat$Time[idx] <= early_window[2],
      stats::runif(m, min = 0.20, max = 0.75),
      stats::runif(m, min = 0.15, max = 0.50)
    )
    shift <- sample(c(-1, 1), size = m, replace = TRUE) * shift_abs
    dat$Time[idx] <- pmin(time_range[2], pmax(time_range[1], dat$Time[idx] + shift))
    already <- dat$Contam_Class[idx] != "none"
    dat$Contam_Class[idx] <- ifelse(already, paste(dat$Contam_Class[idx], "timing_error", sep = "+"), "timing_error")
    dat$Contam_Flag[idx] <- 1L
    dat
  }

  scenario_core <- sub("^SC", "C", scenario)
  out <- switch(
    scenario_core,
    C0_None = out,
    C1_Spike = add_spike(out),
    C2_TimingError = add_timing(out),
    C3_Mixed = add_timing(add_spike(out))
  )

  out$Scenario <- scenario
  out$ScenarioGroup <- if (startsWith(scenario, "SC")) "SparseContamination" else "Contamination"
  rownames(out) <- NULL
  out
}
