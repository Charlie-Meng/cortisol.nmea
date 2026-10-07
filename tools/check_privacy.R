# Privacy guard for the public cortisol.nmea repository.
#
# Fails when tracked (or, with --staged, staged) files include data-like files,
# unusually large files, or identifiers that only occur in private real-data
# material. Run from the repository root:
#   Rscript tools/check_privacy.R          # all tracked files (used in CI)
#   Rscript tools/check_privacy.R --staged # staged files (pre-commit hook)

args <- commandArgs(trailingOnly = TRUE)
staged <- "--staged" %in% args

git_files <- function(staged) {
  cmd <- if (staged) c("diff", "--cached", "--name-only", "--diff-filter=ACMR") else "ls-files"
  out <- system2("git", cmd, stdout = TRUE)
  out[nzchar(out)]
}

files <- git_files(staged)
self <- "tools/check_privacy.R"
allow_path <- "tools/privacy_allowlist.txt"
allowlist <- if (file.exists(allow_path)) {
  x <- trimws(readLines(allow_path, warn = FALSE))
  x[nzchar(x) & !startsWith(x, "#")]
} else character()

blocked_ext <- c("rda", "rds", "rdata", "csv", "tsv", "xlsx", "xls", "sav", "dta",
                 "sas7bdat", "pdf", "pptx", "docx", "feather", "parquet", "sqlite")
max_bytes <- 500 * 1024
# Strings that identify private real-data objects, files, or locations.
private_patterns <- c("ParticipantID", "SpecimenID", "ConvertedConc",
                      "Data_and_Results_preAsy", "Results_postAsy",
                      "real_fitted_parameters", "real_observed_final",
                      "cort_data_time_18", "Dropbox[/\\]Cortisol2021")

problems <- character()
for (f in setdiff(files, c(self, allowlist))) {
  if (!file.exists(f)) next
  ext <- tolower(tools::file_ext(f))
  if (ext %in% blocked_ext) {
    problems <- c(problems, sprintf("%s: data-like file type '.%s'", f, ext))
    next
  }
  size <- file.info(f)$size
  if (is.finite(size) && size > max_bytes) {
    problems <- c(problems, sprintf("%s: file larger than 500 KB (%d bytes)", f, size))
    next
  }
  text <- tryCatch(readLines(f, warn = FALSE, encoding = "UTF-8"), error = function(e) character())
  for (p in private_patterns) {
    if (any(grepl(p, text))) {
      problems <- c(problems, sprintf("%s: contains private identifier pattern '%s'", f, p))
    }
  }
}

if (length(problems)) {
  cat("Privacy check FAILED:\n", paste0("  - ", problems, "\n"), sep = "")
  cat("Remove the file(s), or add a reviewed synthetic file to", allow_path, "\n")
  quit(status = 1)
}
cat(sprintf("Privacy check passed (%d %s files).\n", length(files), if (staged) "staged" else "tracked"))
