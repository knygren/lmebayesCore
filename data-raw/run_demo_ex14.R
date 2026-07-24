devtools::load_all(".", quiet = TRUE)
grDevices::pdf(tempfile(fileext = ".pdf"))
demo("Ex_14_rLMM_estimated_dispersion_estimated_vcov_BigWordClub", package = "lmebayesCore",
     ask = FALSE, echo = FALSE)
grDevices::dev.off()
cat("\nEX_14 DEMO OK\n")
