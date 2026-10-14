#' Write an HTML report of a simulation study
#'
#' Renders tables and figures of a [sim_run()] result: design, failed fits,
#' curve recovery on common subjects, paired differences, an example cohort,
#' outlier detection and subject retention. Requires the rmarkdown package
#' and pandoc.
#'
#' @param runs An `nmea_sim_runs` object.
#' @param file Output HTML file.
#' @param title Report title.
#' @param full_method Name of the workflow method used for paired comparisons,
#'   detection and retention.
#'
#' @return The path of the report, invisibly.
#' @export
sim_report <- function(runs, file = "simulation-report.html",
                       title = "cortisol.nmea simulation report", full_method = "Full") {
  .check_runs(runs)
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package `rmarkdown` is required for sim_report().", call. = FALSE)
  }
  template <- system.file("report", "simulation-report.Rmd", package = "cortisol.nmea")
  work <- tempfile("nmea_report_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  rmd <- file.path(work, "report.Rmd")
  file.copy(template, rmd)
  runs_file <- file.path(work, "runs.rds")
  saveRDS(runs, runs_file)
  out <- rmarkdown::render(rmd, output_file = basename(file), output_dir = work, quiet = TRUE,
                           params = list(runs_file = runs_file, title = title, full_method = full_method),
                           envir = new.env(parent = globalenv()))
  dir.create(dirname(normalizePath(file, mustWork = FALSE)), recursive = TRUE, showWarnings = FALSE)
  file.copy(out, file, overwrite = TRUE)
  invisible(normalizePath(file))
}
