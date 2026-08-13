## Demo: Chapter C07 -- per-group measurement dispersion (rlmerb, bwc)
##
## Same data and formula as Chapter-C04 (bwc full rank), but Block~1 sigma^2_j
## now has a per-school dGamma prior (Ex_13b pattern) instead of fixed pooled
## sigma^2. Block~2 remains dNormal (known variance components).
##
##   demo("Chapter-C07", package = "lmebayesCore")
##
## Legacy reference: demo("Ex_13b_rLMM_estimated_dispersion_known_vcov_BigWordClub_PartVI")
## Vignette Chapter-C07 planned.

.chapter_c07_demo_setup <- function() {
  data_path <- system.file("extdata", "bwc_full_rank.rds",
                           package = "lmebayesCore")
  if (!nzchar(data_path) || !file.exists(data_path)) {
    stop(
      "Cannot find inst/extdata/bwc_full_rank.rds. ",
      "Run data-raw/make_bwc_full_rank.R first.",
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
  stopifnot(all(design$groupef.rank))

  ps <- Prior_Setup_GLMM(
    form,
    data                    = dat,
    pop.pwt                 = 0.01,
    dispformula             = ~school_id,
    group.max_disp_perc     = 0.8,
    group.dispersion.pwt    = 0.1,
    group.alpha_target      = 0.01
  )
  pf <- pfamily_list(ps)
  dl <- dGamma_list(ps, max_disp_perc_measurement = 0.8)

  sizes <- table(dat$school_id)
  mpg   <- attr(dl, "measurement_prior_group")
  sigma2_hat_range <- range(mpg$rate_group / (mpg$shape_group - 1))

  list(
    design          = design,
    ps              = ps,
    pf              = pf,
    dispprior_list  = dl,
    show_groups     = names(sizes)[c(which.max(sizes), which.min(sizes))],
    J               = nlevels(dat$school_id),
    n_obs           = nrow(dat),
    sigma2_hat_range = sigma2_hat_range
  )
}

N_DRAWS <- 1000L
SEED    <- 1L

cat("\n=== Chapter-C07 demo: per-group sigma^2_j (dGamma), rlmerb, n =",
    N_DRAWS, "===\n\n")

S <- .chapter_c07_demo_setup()
cat("Data: bwc full rank, J =", S$J, "schools, n_obs =", S$n_obs, "\n")
cat(sprintf(
  "Per-group sigma^2_hat range [%.3f, %.3f] (dGamma_list from Prior_Setup_GLMM)\n\n",
  S$sigma2_hat_range[1L], S$sigma2_hat_range[2L]
))

set.seed(SEED)
t_fit <- system.time({
  fit <- rlmerb(
    n              = N_DRAWS,
    design         = S$design,
    pfamily_list   = S$pf,
    dispprior_list = S$dispprior_list,
    sim_method     = "TWO_BLOCK_GIBBS",
    progbar        = FALSE,
    verbose        = FALSE,
    print_icm_table = FALSE,
    diag_sweeps    = FALSE
  )
})

cat("Sampler elapsed:", round(t_fit[["elapsed"]], 1), "s\n")
if (!is.null(fit$pilot)) {
  cat("Pilot n:", fit$pilot$n, "\n")
}
cat("m_convergence:", fit$m_convergence, "\n\n")

cat("=== convergence_info ===\n\n")
print(fit$convergence_info)

cat("\n=== summary (two schools: largest and smallest) ===\n\n")
print(summary(fit, groups = S$show_groups))

cat("\n=== plot_var_convergence (main stage, unwhitened) ===\n\n")
plot_var_convergence(fit, whitened = FALSE, stage = "main")

cat("\nDone. See vignette(\"Chapter-C07\") when available.\n")
