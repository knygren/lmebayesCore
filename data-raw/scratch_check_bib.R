## Scratch: confirm every @key cited in the B-series vignettes exists in
## vignettes/REFERENCES.bib.
bib  <- readLines("vignettes/REFERENCES.bib", warn = FALSE)
keys <- sub("^@[A-Za-z]+\\{([^,]+),.*$", "\\1", grep("^@", bib, value = TRUE))
keys <- trimws(sub(",$", "", keys))

cat("--- bib keys ---\n")
print(sort(keys))

files <- list.files("vignettes", pattern = "^Chapter-B.*\\.Rmd$",
                    full.names = TRUE)
ok <- TRUE
for (f in files) {
  txt <- paste(readLines(f, warn = FALSE), collapse = " ")
  m <- regmatches(txt, gregexpr("@[A-Za-z][A-Za-z0-9_]+", txt))[[1]]
  m <- unique(sub("^@", "", m))
  miss <- setdiff(m, keys)
  if (length(miss)) {
    ok <- FALSE
    cat(basename(f), "MISSING:", paste(miss, collapse = ", "), "\n")
  }
}
if (ok) cat("\nAll B-series citations resolve.\n")
