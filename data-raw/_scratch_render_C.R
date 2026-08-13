## Scratch: render the C-series vignettes to data-raw for review.
suppressMessages(devtools::load_all("."))
for (vig in c("vignettes/Chapter-C01.Rmd", "vignettes/Chapter-C02.Rmd",
              "vignettes/Chapter-C03.Rmd", "vignettes/Chapter-C04.Rmd")) {
  nm <- sub("\\.Rmd$", "_preview.html", basename(vig))
  res <- tryCatch({
    rmarkdown::render(vig, output_file = nm, output_dir = "data-raw",
                      quiet = TRUE)
    "OK"
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("%-22s %s\n", basename(vig), res))
}
