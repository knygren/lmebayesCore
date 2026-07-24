devtools::load_all(".", quiet = TRUE)
grDevices::pdf(tempfile(fileext = ".pdf"))
demo("Ex_13_rLMM_estimated_dispersion_known_vcov_BigWordClub", package = "lmebayesCore",
     ask = FALSE, echo = FALSE)
grDevices::dev.off()
cat("\nEX_13 DEMO OK\n")
