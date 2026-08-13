## Append only the words the speller currently flags to inst/WORDLIST,
## preserving the existing curated list and its encoding.
##
## spelling::update_wordlist() is not usable on this machine: it rewrites the
## whole file, drops curated entries, and emits mojibake for UTF-8 sources.

append_missing <- function(pkg) {
  wl_path <- file.path(pkg, "inst", "WORDLIST")

  existing <- readLines(wl_path, encoding = "UTF-8", warn = FALSE)

  found <- spelling::spell_check_package(pkg)
  flagged <- sort(unique(found$word))

  missing <- setdiff(flagged, existing)
  ## Drop mojibake / non-ASCII artefacts rather than baking them in.
  keep <- missing[!grepl("[^A-Za-z'’-]", missing)]
  dropped <- setdiff(missing, keep)

  cat("\n==", basename(pkg), "==\n")
  cat("existing:", length(existing), " flagged:", length(flagged),
      " missing:", length(missing), "\n")
  if (length(dropped)) {
    cat("skipped (non-ASCII / artefact):", paste(dropped, collapse = " "), "\n")
  }
  cat("appending", length(keep), "words:\n")
  print(keep)

  merged <- sort(unique(c(existing, keep)), method = "radix")
  writeLines(merged, wl_path, useBytes = TRUE)
  cat("wrote", length(merged), "words to", wl_path, "\n")

  invisible(keep)
}

append_missing("c:/Rpackages/lmebayes")
append_missing("c:/Rpackages/lmebayesCore")
