## Scratch: diagnose Chapter-C04 vignette / plot failure.
suppressMessages(devtools::load_all("."))
library(lmebayesCore)

dest <- "inst/extdata/Chapter-C04-bwc-twoblock.rds"
cat("artifact exists:", file.exists(dest), "\n")
C04 <- readRDS(dest)
cat("names:", paste(names(C04), collapse = ", "), "\n")
cat("n_chains meta:", C04$meta$n_chains, "\n")
cat("m_convergence:", C04$convergence$m_convergence, "\n")
cat("lambda_star:", C04$convergence$convergence_info$lambda_star, "\n")

sh <- C04$sweep_history
cat("sweep_history stage:", sh$stage, " n_sweeps:", sh$n_sweeps, "\n")
cat("table rows:", nrow(sh$table), "\n")

cat("\nexact ref test:\n")
ref <- try(lmerb_posterior_covariance(C04$design, C04$measurement_prior_list),
           silent = TRUE)
if (inherits(ref, "try-error")) {
  cat("FAILED:", conditionMessage(attr(ref, "condition")), "\n")
} else {
  cat("OK dim", paste(dim(ref), collapse = "x"), "\n")
}

cat("\nplot_var_convergence test:\n")
res <- try(plot_var_convergence(
  C04$sweep_history,
  design                 = C04$design,
  measurement_prior_list   = C04$measurement_prior_list,
  whitened               = FALSE,
  n_chains               = C04$meta$n_chains
), silent = TRUE)
if (inherits(res, "try-error")) {
  cat("PLOT FAILED:", conditionMessage(attr(res, "condition")), "\n")
  print(attr(res, "condition"))
} else {
  cat("plot OK\n")
}

cat("\nsummary object class:", class(C04$summary)[1], "\n")
