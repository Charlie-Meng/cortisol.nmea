# Installs a local git pre-commit hook that runs the privacy check on staged files.
# Run once from the repository root with the R you develop with:
#   Rscript tools/install_hooks.R
rscript <- normalizePath(file.path(R.home("bin"), "Rscript"), winslash = "/", mustWork = FALSE)
hook <- file.path(".git", "hooks", "pre-commit")
writeLines(c("#!/bin/sh", sprintf("\"%s\" tools/check_privacy.R --staged || exit 1", rscript)), hook)
Sys.chmod(hook, "755")
cat("Installed", hook, "using", rscript, "\n")
