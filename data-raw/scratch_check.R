## Scratch: run R CMD check (vignettes included) and print a compact summary.
## Usage: Rscript data-raw/scratch_check.R [path-to-package]
args <- commandArgs(trailingOnly = TRUE)
pkg  <- if (length(args)) args[[1L]] else "."

r <- rcmdcheck::rcmdcheck(
  path     = pkg,
  args     = c("--no-manual", "--as-cran"),
  error_on = "never"
)

cat("\n\n===================== SUMMARY:", pkg, "=====================\n")
cat("--- ERRORS (", length(r$errors), ") ---\n", sep = "")
if (length(r$errors))   cat(r$errors,   sep = "\n---\n")
cat("--- WARNINGS (", length(r$warnings), ") ---\n", sep = "")
if (length(r$warnings)) cat(r$warnings, sep = "\n---\n")
cat("--- NOTES (", length(r$notes), ") ---\n", sep = "")
if (length(r$notes))    cat(r$notes,    sep = "\n---\n")
cat("=============================================================\n")
