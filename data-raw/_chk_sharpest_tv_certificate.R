## Scratch validation for gamma_beta_tv_certificate() pipeline.
## Not part of the package test suite.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required to load lmebayesCore for this script.")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Suggested package 'bayesrules' is required for the Gaussian fixture.")
}

devtools::load_all("c:/Rpackages/lmebayesCore", quiet = TRUE)

cat("=== Sharpest TV certificate: Gaussian smoke test ===\n")

data(big_word_club, package = "bayesrules")
dat_g <- big_word_club
dat_g$school_id <- factor(dat_g$school_id)
dat_g <- subset(
  dat_g,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat_g[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

design0 <- model_setup(form_lmer, data = dat_g)
dat_g <- subset(dat_g, school_id %in% names(design0$groupef.rank)[design0$groupef.rank])
dat_g$school_id <- droplevels(dat_g$school_id)

design <- model_setup(form_lmer, data = dat_g)
ps <- Prior_Setup_GLMM(form_lmer, data = dat_g, pop.pwt = 0.01)
pf <- pfamily_list(ps)
dispprior_list <- list(dispersion = ps$group.dispersion)

cert <- gamma_beta_tv_certificate(
  design = design,
  pfamily_list = pf,
  family = gaussian(),
  delta_2 = 0.01,
  dispprior_list = dispprior_list,
  k = 50L,
  kappa_method = "none",
  estep = "exact"
)

stopifnot(inherits(cert, "gamma_beta_tv_certificate"))
stopifnot(isTRUE(cert$certified$gamma_em))
stopifnot(isTRUE(cert$certified$epsilon))
stopifnot(identical(cert$certified$sharpest_display, "limit"))
stopifnot(cert$rosenthal$bound > 0)
stopifnot(cert$rosenthal$bound < 1)
stopifnot(cert$full_pi_gamma$full_bound >= cert$full_pi_gamma$inner_bound)
stopifnot(abs(cert$full_pi_gamma$full_bound - cert$full_pi_gamma$inner_bound - cert$delta_2) < 1e-12)

cat("  inner bound at k=50:", signif(cert$rosenthal$bound, 4), "\n")
cat("  full pi_gamma bound:", signif(cert$full_pi_gamma$full_bound, 4), "\n")
cat("  kappa_max^LB:", signif(cert$floor_spectrum$kappa_max_lb, 4), "\n")
print(cert)

cat("\n=== Floor spectrum + Rosenthal sanity (from certificate) ===\n")

safe <- cert$beta_set
stopifnot(inherits(safe, "beta_marginal_safe_set"))
stopifnot(identical(safe$level$scheme, "r_gauss_joint"))
stopifnot(safe$level$r_gauss > 0)

spec1 <- floor_coupling_spectrum(cert$mode, safe)

cat("\n=== Rosenthal alpha optimization sanity ===\n")

ros_opt <- rosenthal_tv_bound(
  k = 100L,
  eps = cert$epsilon$eps_star,
  spectrum = spec1,
  display_mode = "sharp"
)
ros_fix <- rosenthal_tv_bound(
  k = 100L,
  eps = cert$epsilon$eps_star,
  spectrum = spec1,
  alpha = ros_opt$alpha,
  display_mode = "sharp"
)
stopifnot(abs(ros_opt$bound - ros_fix$bound) < 1e-12)
stopifnot(ros_opt$alpha > spec1$lambda_lb && ros_opt$alpha < 1)

cat("  optimized alpha:", signif(ros_opt$alpha, 6), "\n")
cat("  bound:", signif(ros_opt$bound, 6), "\n")

cat("\nAll sharpest TV certificate checks passed.\n")
