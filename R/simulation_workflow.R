#' Run a synthetic cortisol simulation study
#'
#' High-level public workflow that generates latent NMEA curves, observed
#' synthetic measurements, missing/sparse scenarios, and contamination
#' scenarios. The output is designed for package examples and benchmark
#' pipelines and contains no real participant data.
#'
#' @param n Number of synthetic subjects.
#' @param time_grid Latent curve evaluation grid.
#' @param theta Mean parameter vector for [simulate_nmea_parameters()].
#' @param covariance Covariance matrix for [simulate_nmea_parameters()].
#' @param alpha Fixed alpha value for synthetic subjects.
#' @param time_templates Time templates passed to [simulate_observed_data()].
#' @param sigma Scalar or function passed to [simulate_observed_data()].
#' @param missing_scenarios Missing/sparse scenarios to generate.
#' @param contamination_scenarios Contamination scenarios to generate.
#' @param seed Optional random seed.
#'
#' @return A list with synthetic parameters, latent curves, observed data,
#'   scenario data, and configuration metadata.
#' @export
run_simulation_study <- function(n = 50,
                                 time_grid = seq(0, 18, by = 0.25),
                                 theta = c(mu = 1.0, s = 1.2, c1 = 35, c0 = 2),
                                 covariance = diag(c(0.25, 0.12, 5, 0.5)^2),
                                 alpha = 1.15,
                                 time_templates = list(
                                   standard = c(0, 0.5, 1, 3, 6, 9, 12, 15),
                                   sparse = c(0, 0.75, 3, 8, 14)
                                 ),
                                 sigma = function(time) 0.7 + 0.15 * sqrt(pmax(time, 0)),
                                 missing_scenarios = c(
                                   "S0_CleanObserved",
                                   "S1_MCAR15",
                                   "S2_MorningMissing",
                                   "S3_Sparse34",
                                   "S4_SparseMorningMissing"
                                 ),
                                 contamination_scenarios = c(
                                   "C0_None",
                                   "C1_Spike",
                                   "C2_TimingError",
                                   "C3_Mixed",
                                   "SC0_None",
                                   "SC1_Spike",
                                   "SC2_TimingError",
                                   "SC3_Mixed"
                                 ),
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  seed_base <- if (is.null(seed)) sample.int(.Machine$integer.max, 1) else seed

  psi <- simulate_nmea_parameters(
    n = n,
    theta = theta,
    covariance = covariance,
    alpha = alpha,
    seed = seed_base
  )
  latent <- simulate_latent_curves(psi, time_grid = time_grid)
  observed <- simulate_observed_data(
    psi = psi,
    time_templates = time_templates,
    sigma = sigma,
    seed = seed_base + 1L
  )

  scenarios <- list()
  for (i in seq_along(missing_scenarios)) {
    sc <- missing_scenarios[i]
    scenarios[[sc]] <- apply_missing_scenario(
      observed,
      scenario = sc,
      seed = seed_base + 100L + i
    )
  }

  base_clean <- scenarios[["S0_CleanObserved"]]
  if (is.null(base_clean)) {
    base_clean <- apply_missing_scenario(observed, "S0_CleanObserved", seed = seed_base + 10L)
  }
  base_sparse <- scenarios[["S4_SparseMorningMissing"]]
  if (is.null(base_sparse)) {
    base_sparse <- apply_missing_scenario(observed, "S4_SparseMorningMissing", seed = seed_base + 11L)
  }

  for (i in seq_along(contamination_scenarios)) {
    sc <- contamination_scenarios[i]
    base <- if (startsWith(sc, "SC")) base_sparse else base_clean
    scenarios[[sc]] <- apply_contamination_scenario(
      base,
      scenario = sc,
      seed = seed_base + 200L + i
    )
  }

  structure(
    list(
      psi = psi,
      latent = latent,
      observed = observed,
      scenarios = scenarios,
      config = list(
        n = n,
        time_grid = time_grid,
        alpha = alpha,
        missing_scenarios = missing_scenarios,
        contamination_scenarios = contamination_scenarios,
        seed = seed
      )
    ),
    class = "nmea_simulation"
  )
}

#' Run method comparisons on a synthetic simulation object
#'
#' Compares one or more methods against the known latent truth in a
#' `run_simulation_study()` object. Fast truth-based methods are available for
#' vignettes and tests; fitted methods are available when optional dependencies
#' are installed.
#'
#' @param simulation Object returned by [run_simulation_study()].
#' @param methods Character vector of methods. Supported values are
#'   `OracleTruth`, `SymmetricTruth`, `Baseline_NMEA`, `NMEA_AlphaOnly`, and
#'   `GAM`.
#' @param scenarios Scenario names to evaluate.
#' @param alpha_grid Alpha grid for `NMEA_AlphaOnly`.
#' @param fit_control Optional `saemix` control list for fitted NMEA methods.
#' @param seed Optional seed used by fitted methods.
#' @param keep_predictions Keep prediction matrices in the returned object.
#'
#' @return A list with method-level summary rows, subject-level metrics, and
#'   optionally prediction matrices.
#' @export
run_method_comparison <- function(simulation,
                                  methods = c("OracleTruth", "SymmetricTruth"),
                                  scenarios = names(simulation$scenarios),
                                  alpha_grid = seq(0.8, 1.5, by = 0.1),
                                  fit_control = list(),
                                  seed = 632545,
                                  keep_predictions = FALSE) {
  if (!inherits(simulation, "nmea_simulation")) {
    stop("`simulation` must be returned by `run_simulation_study()`.", call. = FALSE)
  }

  methods <- match.arg(
    methods,
    choices = c("OracleTruth", "SymmetricTruth", "Baseline_NMEA", "NMEA_AlphaOnly", "GAM"),
    several.ok = TRUE
  )

  time_grid <- simulation$latent$time_grid
  truth_all <- simulation$latent$curve_mat
  summary_rows <- list()
  subject_rows <- list()
  prediction_rows <- list()

  for (scenario in scenarios) {
    dat <- simulation$scenarios[[scenario]]
    if (is.null(dat)) next
    subject_ids <- intersect(unique(as.character(dat$Subject)), rownames(truth_all))
    truth <- truth_all[subject_ids, , drop = FALSE]

    for (method in methods) {
      pred_obj <- .run_one_method(
        method = method,
        dat = dat,
        psi = simulation$psi,
        subject_ids = subject_ids,
        time_grid = time_grid,
        alpha_grid = alpha_grid,
        fit_control = fit_control,
        seed = seed
      )

      metric_df <- if (!is.null(pred_obj$pred)) {
        curve_benchmark_metrics(pred_obj$pred, truth, time_grid)
      } else {
        data.frame(Subject = subject_ids, Curve_RMSE = NA_real_, AUC_Error = NA_real_)
      }

      sum_df <- summarize_curve_metrics(metric_df)
      sum_df$Scenario <- scenario
      sum_df$ScenarioGroup <- unique(dat$ScenarioGroup)[1]
      sum_df$Method <- method
      sum_df$Status <- pred_obj$status
      sum_df$Failure_Reason <- pred_obj$failure_reason
      summary_rows[[paste(scenario, method, sep = "::")]] <- sum_df

      metric_df$Scenario <- scenario
      metric_df$ScenarioGroup <- unique(dat$ScenarioGroup)[1]
      metric_df$Method <- method
      subject_rows[[paste(scenario, method, sep = "::")]] <- metric_df

      if (isTRUE(keep_predictions)) {
        prediction_rows[[paste(scenario, method, sep = "::")]] <- pred_obj$pred
      }
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  subject_df <- do.call(rbind, subject_rows)
  rownames(subject_df) <- NULL

  structure(
    list(
      summary = summary_df,
      subject_metrics = subject_df,
      predictions = if (isTRUE(keep_predictions)) prediction_rows else NULL
    ),
    class = "nmea_method_comparison"
  )
}

.run_one_method <- function(method,
                            dat,
                            psi,
                            subject_ids,
                            time_grid,
                            alpha_grid,
                            fit_control,
                            seed) {
  ok <- function(pred) list(pred = pred, status = "ok", failure_reason = NA_character_)
  fail <- function(e) list(pred = NULL, status = "failed", failure_reason = conditionMessage(e))

  tryCatch({
    if (identical(method, "OracleTruth")) {
      return(ok(model_eval(time_grid, psi[subject_ids, , drop = FALSE])))
    }

    if (identical(method, "SymmetricTruth")) {
      psi2 <- psi[subject_ids, , drop = FALSE]
      psi2[, "alpha"] <- 1
      return(ok(model_eval(time_grid, psi2)))
    }

    if (identical(method, "Baseline_NMEA")) {
      fit <- fit_baseline_nmea(
        dat,
        id = "Subject",
        time = "Time",
        response = "Conc_Obs",
        seed = seed,
        control = fit_control
      )
      return(ok(predict_nmea(fit, time_grid)))
    }

    if (identical(method, "NMEA_AlphaOnly")) {
      fit <- fit_alpha_only_nmea(
        dat,
        alpha_grid = alpha_grid,
        id = "Subject",
        time = "Time",
        response = "Conc_Obs",
        seed = seed,
        control = fit_control
      )
      return(ok(predict_nmea(fit, time_grid)))
    }

    if (identical(method, "GAM")) {
      fit <- fit_gam_baseline(dat, id = "Subject", time = "Time", response = "Conc_Obs")
      return(ok(predict_gam_baseline(fit, subject_ids = subject_ids, time_grid = time_grid)))
    }

    stop("Unsupported method: ", method, call. = FALSE)
  }, error = fail)
}
