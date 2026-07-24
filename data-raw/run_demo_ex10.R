devtools::load_all(".", quiet = TRUE)
grDevices::pdf(tempfile(fileext = ".pdf"))
demo("Ex_10_rLMM_known_dispersion_known_vcov_BigWordClub", package = "lmebayesCore",
     ask = FALSE, echo = FALSE)
grDevices::dev.off()
cat("\nEX_10 DEMO OK\n")
