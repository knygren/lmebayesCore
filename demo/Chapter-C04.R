## Demo: Chapter C04 -- two-block Gibbs sampler (rlmerb TWO_BLOCK_GIBBS)
##
## Interactive counterpart to vignette("Chapter-C04"): same bwc model and priors
## as Chapter-C02, n = 1000 draws (vignette stores n = 10000).
## progbar = FALSE: the sweep progress bar floods console capture on long runs.
##
##   demo("Chapter-C04", package = "lmebayesCore")
##
## Manual test: tests/manual/test-manual-Chapter-C04.R
## Artifact:    data-raw/make_Chapter-C04.R -> inst/extdata/Chapter-C04-bwc-twoblock.rds

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

cat("\n=== Chapter-C04 demo: two-block Gibbs (TWO_BLOCK_GIBBS), n =",
    N_DRAWS, "===\n\n")

S <- .chapter_c_demo_setup()

set.seed(SEED)
t_fit <- system.time({
  fit <- rlmerb(
    n              = N_DRAWS,
    design         = S$design,
    pfamily_list   = S$pf,
    dispprior_list = list(dispersion = S$ps$group.dispersion),
    sim_method     = "TWO_BLOCK_GIBBS",
    progbar        = FALSE,
    verbose        = FALSE,
    print_icm_table = FALSE,
    diag_sweeps    = FALSE
  )
})

cat("\nSampler elapsed:", round(t_fit[["elapsed"]], 1), "s\n")
cat("m_convergence:", fit$m_convergence, "\n\n")

cat("=== convergence_info ===\n\n")
print(fit$convergence_info)

cat("\n=== summary (two schools: largest and smallest) ===\n\n")
print(summary(fit, groups = S$show_groups))

cat("\nDone. See vignette(\"Chapter-C04\") for sweep history and theory.\n")
