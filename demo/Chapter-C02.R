## Demo: Chapter C02 -- exact Gaussian iid sampler (rlmerb DEFAULT)
##
## Interactive counterpart to vignette("Chapter-C02"): same bwc model and priors
## as Chapter-C01, n = 1000 draws (vignette stores n = 10000).
##
##   demo("Chapter-C02", package = "lmebayesCore")
##
## Manual test: tests/manual/test-manual-Chapter-C02.R (n = 10000)
## Artifact:    data-raw/make_Chapter-C02.R -> inst/extdata/Chapter-C02-bwc-iid.rds

.chapter_c_demo_setup <- function() {
  data_path <- system.file("extdata", "bwc_full_rank.rds",
                           package = "lmebayesCore")
  if (!nzchar(data_path) || !file.exists(data_path)) {
    stop(
      "Cannot find inst/extdata/bwc_full_rank.rds. ",
      "Reinstall lmebayesCore or run data-raw/make_bwc_full_rank.R.",
      call. = FALSE
    )
  }
  dat  <- readRDS(data_path)
  form <- score_ppvt ~
    private_school + title1 + free_reduced_lunch +
    distracted_ppvt + distracted_a1 +
    free_reduced_lunch:distracted_a1 +
    (1 + distracted_ppvt + distracted_a1 || school_id)
  design <- model_setup(form, data = dat)
  ps     <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)
  pf     <- pfamily_list(ps)
  sizes  <- table(dat$school_id)
  list(
    design      = design,
    ps          = ps,
    pf          = pf,
    show_groups = names(sizes)[c(which.max(sizes), which.min(sizes))]
  )
}

N_DRAWS <- 1000L
SEED    <- 1L

cat("\n=== Chapter-C02 demo: exact iid (sim_method = DEFAULT), n =",
    N_DRAWS, "===\n\n")

S <- .chapter_c_demo_setup()
cat("Data: bayesrules big_word_club (full rank), J =",
    nlevels(S$design$group), "schools\n\n")

set.seed(SEED)
t_fit <- system.time({
  fit <- rlmerb(
    n              = N_DRAWS,
    design         = S$design,
    pfamily_list   = S$pf,
    dispprior_list = list(dispersion = S$ps$group.dispersion),
    sim_method     = "DEFAULT",
    progbar        = FALSE,
    verbose        = FALSE,
    print_icm_table = FALSE
  )
})

cat("Sampler elapsed:", round(t_fit[["elapsed"]], 1), "s\n\n")

cat("=== convergence_info ===\n\n")
print(fit$convergence_info)

cat("\n=== summary (two schools: largest and smallest) ===\n\n")
print(summary(fit, groups = S$show_groups))

cat("\nDone. See vignette(\"Chapter-C02\") for the full n = 10000 analysis.\n")
