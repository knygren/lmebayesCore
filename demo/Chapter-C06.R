## Demo: Chapter C06 -- binomial GLMM (rglmerb, book_banning)
##
## Same model as inst/examples/Ex_rglmerb.R; n = 1000 (help example uses n = 10).
## progbar = FALSE: quiet console for scripted runs; enable later for interactive demo().
##
##   demo("Chapter-C06", package = "lmebayesCore")
##
## Vignette Chapter-C06 planned. Manual test: tests/manual/ (TBD).

.chapter_c06_demo_setup <- function() {
  if (!requireNamespace("bayesrules", quietly = TRUE)) {
    stop("This demo requires the 'bayesrules' package.", call. = FALSE)
  }
  data(book_banning, package = "bayesrules")
  dat <- book_banning[, c("state", "removed", "violent")]
  dat <- dat[stats::complete.cases(dat), ]
  dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
  dat$violent_i <- as.integer(
    dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
  )
  dat$state <- factor(dat$state)
  keep <- names(sort(table(dat$state), decreasing = TRUE))[seq_len(8L)]
  dat <- droplevels(subset(dat, state %in% keep))

  form <- removed_i ~ violent_i + (1 + violent_i || state)
  design <- model_setup(form, data = dat, family = binomial())
  ps <- Prior_Setup_GLMM(form, data = dat, family = binomial(), pop.pwt = 0.01)

  list(
    design = design,
    ps     = ps,
    pf     = pfamily_list(ps),
    J      = nlevels(dat$state),
    n_obs  = nrow(dat)
  )
}

N_DRAWS <- 1000L
SEED    <- 1L

cat("\n=== Chapter-C06 demo: binomial GLMM (rglmerb), n =", N_DRAWS, "===\n\n")

S <- .chapter_c06_demo_setup()
cat("Data: book_banning (8 largest states), J =", S$J,
    "states, n_obs =", S$n_obs, "\n\n")

set.seed(SEED)
t_fit <- system.time({
  fit <- rglmerb(
    n              = N_DRAWS,
    design         = S$design,
    pfamily_list   = S$pf,
    family         = binomial(),
    dispprior_list = NULL,
    progbar        = TRUE,
    verbose        = TRUE
  )
})

cat("Sampler elapsed:", round(t_fit[["elapsed"]], 1), "s\n")
if (!is.null(fit$pilot)) {
  cat("Pilot n:", fit$pilot$n, "\n")
}
cat("m_convergence:", fit$m_convergence, "\n\n")

cat("=== convergence_info ===\n\n")
print(fit$convergence_info)

cat("\n=== summary ===\n\n")
print(summary(fit))

cat("\n=== print_groupef (draws 1:3, groups 1:3) ===\n\n")
print_groupef(fit, draws = 1:3, groups = 1:3)

cat("\nDone. Vignette Chapter-C06 planned.\n")
