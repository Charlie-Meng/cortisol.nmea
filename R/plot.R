#' Plot simulation results
#'
#' ggplot2 figures reproducing the panels of the paper's simulation figures.
#'
#' * `plot_cohort_curves()`: individual fitted curves (thin) and each method's
#'   mean curve (thick) with the mean true curve (dashed), one panel per
#'   scenario for one base cohort. Mean curves use each method's own available
#'   subjects.
#' * `plot_subject()`: one subject's fitted and true curves with its
#'   measurements classified by the workflow's outlier step: TP (contaminated,
#'   flagged), FP (clean, flagged), FN (contaminated, missed), TN (clean, kept).
#' * `plot_recovery()`: cohort-level curve RMSE or curve-AUC error on common
#'   subjects (small dots), with means and 95% Monte Carlo intervals.
#' * `plot_detection()`: FPR and FNR of a workflow method by trimester.
#' * `plot_retention()`: stacked percentages of retained and excluded subjects.
#'
#' @param runs An `nmea_sim_runs` object.
#' @param trimester,seed Identify the base cohort (default: first in design).
#' @param scenarios Scenarios to show.
#' @param scenario Scenario of the subject.
#' @param subject Subject identifier.
#' @param methods Methods to show (default: all).
#' @param flag_method Workflow method whose flags classify the measurements.
#' @param y_max Upper limit of the y axis for cohort curves (`NULL`: automatic).
#' @param metric `"rmse"` or `"auc_error"`.
#' @param method Workflow method for detection or retention.
#'
#' @return A ggplot object.
#' @name plot_simulation
NULL

.method_colours <- function(methods) {
  base <- c(Direct = "#23866C", Full = "#277DA8", GAMM = "#C05C35")
  extra <- setdiff(methods, names(base))
  if (length(extra)) {
    pal <- c("#8267AB", "#B8860B", "#4B5563", "#A23B72", "#3B8EA5", "#6B8E23")
    base <- c(base, stats::setNames(rep_len(pal, length(extra)), extra))
  }
  c(base[methods], Truth = "#111111")
}
.point_colours <- c(TN = "#4B5563", TP = "#D55E00", FP = "#CC79A7", FN = "#E69F00")
.point_shapes <- c(TN = 1, TP = 17, FP = 4, FN = 15)

.theme_nmea <- function() {
  ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.28, colour = "#374151"),
      strip.background = ggplot2::element_rect(fill = "#F1F3F5", colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 9.5),
      legend.position = "bottom"
    )
}

.find_case <- function(runs, trimester, seed, scenario) {
  id <- sprintf("%s_%s_%s", trimester, seed, scenario)
  cs <- runs$cases[[id]]
  if (is.null(cs)) stop("No case ", id, " in `runs`.", call. = FALSE)
  cs
}

.curve_long <- function(mat, times, method) {
  data.frame(subject = rep(rownames(mat), each = length(times)), time = rep(times, nrow(mat)),
             value = as.vector(t(mat)), method = method, stringsAsFactors = FALSE)
}

#' @rdname plot_simulation
#' @export
plot_cohort_curves <- function(runs, trimester = runs$design$trimesters[1],
                               seed = runs$design$seeds[[trimester]][1],
                               scenarios = names(runs$design$scenarios), methods = NULL,
                               y_max = 20) {
  .check_runs(runs)
  methods <- .check_methods(runs, methods)
  times <- runs$design$times
  ind <- list(); avg <- list()
  for (sc in scenarios) {
    cs <- .find_case(runs, trimester, seed, sc)
    for (m in methods) {
      cur <- cs$methods[[m]]$curves
      ok <- apply(cur, 1, function(x) all(is.finite(x)))
      if (any(ok)) {
        ind[[length(ind) + 1]] <- cbind(.curve_long(cur[ok, , drop = FALSE], times, m), scenario = sc)
        avg[[length(avg) + 1]] <- data.frame(time = times, value = colMeans(cur[ok, , drop = FALSE]),
                                             method = m, scenario = sc, n = sum(ok))
      }
    }
    avg[[length(avg) + 1]] <- data.frame(time = times, value = colMeans(cs$truth_curves),
                                         method = "Truth", scenario = sc, n = nrow(cs$truth_curves))
  }
  ind <- do.call(rbind, ind); avg <- do.call(rbind, avg)
  lev <- c(methods, "Truth")
  ind$method <- factor(ind$method, lev); avg$method <- factor(avg$method, lev)
  counts <- stats::aggregate(n ~ scenario + method, avg[avg$method != "Truth", ], function(x) x[1])
  lab <- vapply(scenarios, function(sc) {
    z <- counts[counts$scenario == sc, ]
    sprintf("%s (n: %s)", sc, paste(z$n[order(match(z$method, methods))], collapse = "/"))
  }, character(1))
  ind$scenario <- factor(ind$scenario, scenarios, lab); avg$scenario <- factor(avg$scenario, scenarios, lab)
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = ind, ggplot2::aes(.data$time, .data$value, colour = .data$method,
                                                group = interaction(.data$method, .data$subject)),
                       alpha = 0.12, linewidth = 0.15) +
    ggplot2::geom_line(data = avg, ggplot2::aes(.data$time, .data$value, colour = .data$method,
                                                linetype = .data$method), linewidth = 0.8) +
    ggplot2::facet_wrap(~scenario, nrow = 1) +
    ggplot2::scale_colour_manual(values = .method_colours(methods), name = NULL) +
    ggplot2::scale_linetype_manual(values = c(stats::setNames(rep("solid", length(methods)), methods),
                                              Truth = "22"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = c(0, 6, 12, 18)) +
    ggplot2::labs(x = "Hours since waking", y = "Cortisol (nmol/L)",
                  title = sprintf("%s, seed %s: individual and mean curves", trimester, seed)) +
    .theme_nmea()
  if (!is.null(y_max)) p <- p + ggplot2::coord_cartesian(ylim = c(0, y_max))
  p
}

#' @rdname plot_simulation
#' @export
plot_subject <- function(runs, subject, trimester = runs$design$trimesters[1],
                         seed = runs$design$seeds[[trimester]][1], scenario = "D2",
                         methods = NULL, flag_method = "Full") {
  .check_runs(runs)
  methods <- .check_methods(runs, methods)
  cs <- .find_case(runs, trimester, seed, scenario)
  if (!subject %in% rownames(cs$truth_curves)) stop("Unknown subject ", subject, call. = FALSE)
  times <- runs$design$times
  cur <- do.call(rbind, lapply(methods, function(m) {
    .curve_long(cs$methods[[m]]$curves[subject, , drop = FALSE], times, m)
  }))
  cur <- rbind(cur, .curve_long(cs$truth_curves[subject, , drop = FALSE], times, "Truth"))
  cur$method <- factor(cur$method, c(methods, "Truth"))
  pts <- cs$points[cs$points$subject == subject, ]
  pts$time_observed <- cs$observed_time[match(pts$point_id, cs$points$point_id)]
  fl <- cs$methods[[flag_method]]$points
  flagged <- if (is.null(fl)) rep(FALSE, nrow(pts)) else fl$flagged[match(pts$point_id, fl$point_id)]
  pts$class <- factor(ifelse(pts$contaminated, ifelse(flagged, "TP", "FN"), ifelse(flagged, "FP", "TN")),
                      names(.point_shapes))
  outcome <- cs$methods[[flag_method]]$outcome[[subject]]
  line_keys <- c(methods, "Truth")
  # One colour scale for lines and points; the colour legend lists only the lines,
  # and the shape legend carries the point colours.
  colours <- c(.method_colours(methods), .point_colours)
  ggplot2::ggplot() +
    ggplot2::geom_line(data = cur, ggplot2::aes(.data$time, .data$value, colour = .data$method,
                                                linetype = .data$method), linewidth = 0.7, na.rm = TRUE) +
    ggplot2::geom_point(data = pts, ggplot2::aes(.data$time_observed, .data$y, shape = .data$class,
                                                 colour = .data$class),
                        size = 2, stroke = 0.7, show.legend = c(colour = FALSE, shape = TRUE)) +
    ggplot2::scale_colour_manual(values = colours, breaks = line_keys, name = NULL) +
    ggplot2::scale_linetype_manual(values = c(stats::setNames(rep("solid", length(methods)), methods),
                                              Truth = "22"), breaks = line_keys, name = NULL) +
    ggplot2::scale_shape_manual(values = .point_shapes, drop = FALSE, name = NULL,
                                labels = c(TN = "TN: clean, kept", TP = "TP: contaminated, flagged",
                                           FP = "FP: clean, flagged", FN = "FN: contaminated, missed")) +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = unname(.point_colours)),
                                                  order = 2),
                    colour = ggplot2::guide_legend(order = 1), linetype = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_x_continuous(breaks = c(0, 6, 12, 18), limits = c(0, 18)) +
    ggplot2::labs(x = "Hours since waking", y = "Cortisol (nmol/L)",
                  title = sprintf("%s | %s (%s: %s)", subject, scenario, flag_method, outcome)) +
    .theme_nmea() +
    ggplot2::theme(legend.box = "vertical", legend.spacing.y = ggplot2::unit(0, "mm"))
}

#' @rdname plot_simulation
#' @export
plot_recovery <- function(runs, metric = c("rmse", "auc_error"), methods = NULL) {
  metric <- match.arg(metric)
  methods <- .check_methods(runs, methods)
  s <- eval_summary(runs, methods)
  coh <- s$cohorts[s$cohorts$valid, ]
  sm <- s$summary[s$summary$measure == metric & s$summary$trimester != "All", ]
  scen <- names(runs$design$scenarios)
  offset <- function(m) (match(m, methods) - (length(methods) + 1) / 2) * 0.22
  coh$x <- match(coh$scenario, scen) + offset(coh$method)
  sm$x <- match(sm$scenario, scen) + offset(sm$method)
  coh$value <- coh[[metric]]
  coh$method <- factor(coh$method, methods); sm$method <- factor(sm$method, methods)
  ylab <- if (metric == "rmse") "Curve RMSE (nmol/L)" else "Curve-AUC error (nmol h/L)"
  ggplot2::ggplot() +
    ggplot2::geom_point(data = coh, ggplot2::aes(.data$x, .data$value, colour = .data$method),
                        size = 0.9, alpha = 0.3, na.rm = TRUE) +
    ggplot2::geom_errorbar(data = sm, ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper,
                                                   colour = .data$method), width = 0.1, na.rm = TRUE) +
    ggplot2::geom_point(data = sm, ggplot2::aes(.data$x, .data$mean, colour = .data$method), size = 2,
                        na.rm = TRUE) +
    ggplot2::facet_wrap(~trimester, nrow = 1) +
    ggplot2::scale_colour_manual(values = .method_colours(methods), name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq_along(scen), labels = scen) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = NULL, y = ylab,
                  title = "Curve recovery on common subjects (dots: cohorts; bars: 95% MC interval)") +
    .theme_nmea()
}

#' @rdname plot_simulation
#' @export
plot_detection <- function(runs, method = "Full", scenario = "D2") {
  d <- eval_detection(runs, method)$cohorts
  d <- d[d$scenario == scenario, ]
  long <- rbind(data.frame(trimester = d$trimester, measure = "FPR", value = 100 * d$FPR),
                data.frame(trimester = d$trimester, measure = "FNR", value = 100 * d$FNR))
  long$measure <- factor(long$measure, c("FPR", "FNR"))
  ggplot2::ggplot(long, ggplot2::aes(.data$trimester, .data$value, colour = .data$measure,
                                     shape = .data$measure)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.4), alpha = 0.4, na.rm = TRUE) +
    ggplot2::stat_summary(fun = mean, geom = "point", size = 3,
                          position = ggplot2::position_dodge(width = 0.4), na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = c(FPR = "#CC79A7", FNR = "#E69F00"),
                                 labels = c("FPR: clean points flagged", "FNR: contamination missed"), name = NULL) +
    ggplot2::scale_shape_manual(values = c(FPR = 4, FNR = 15),
                                labels = c("FPR: clean points flagged", "FNR: contamination missed"), name = NULL) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = NULL, y = "Point-level error (%)",
                  title = sprintf("Detection errors of %s, %s", method, scenario)) +
    .theme_nmea()
}

#' @rdname plot_simulation
#' @export
plot_retention <- function(runs, method = "Full", scenario = "D2") {
  s <- eval_retention(runs, method)$summary
  s <- s[s$scenario == scenario & s$trimester != "All", ]
  labels <- c(retained = "Retained", excluded_min_obs = "Excluded: too few points",
              excluded_fvu = "Excluded: FVU", excluded_c1 = "Excluded: c1 <= 0",
              not_fitted = "No curve")
  s$measure <- factor(s$measure, names(labels), labels)
  s <- s[is.finite(s$mean) & s$mean > 0, ]
  ggplot2::ggplot(s, ggplot2::aes(.data$trimester, .data$mean, fill = .data$measure)) +
    ggplot2::geom_col(width = 0.55, position = ggplot2::position_stack(reverse = TRUE),
                      colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = stats::setNames(c("#277DA8", "#9AA7B3", "#D3D9DE", "#E7EBEF", "#F3F4F6"),
                                                        labels), name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 100.5), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Subjects (%)",
                  title = sprintf("Subject retention of %s, %s", method, scenario)) +
    .theme_nmea()
}
