# Scratch checks for restricted-chain epsilon (Chapter C05 / MVN route).
# Not part of the package test suite.
#
# Gaussian parity (fast): data-raw/_ex_c05_gaussian_parity.R (big_word_club).
# This file: MVN arithmetic + big_word_club Gaussian certificate smoke test.
# Binomial MC checks are slow; run separately when needed.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required to load lmebayesCore for this script.")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop(
    "Suggested package 'bayesrules' is required ",
    "(same datasets as inst/examples/Ex_rglmerb.R and Ex_rLMM_reg.R)."
  )
}
devtools::load_all(".", quiet = TRUE)

cat("=== Check A: legacy mvn_calibrate + equal-kappa spectrum parity ===\n")

q <- 3L
delta <- 0.01
tol <- 0.05

cal <- mvn_calibrate(q, delta)
stopifnot(abs(cal$d - 0.5 * cal$r2) < 1e-12)
stopifnot(abs(cal$set_cut - exp(-cal$d)) < 1e-15)

spec_equal <- list(
  kappa = rep(0.5, q),
  weights = rep(1, q),
  kappa_max = 0.5,
  q = q,
  method = "symmetrized_coupling"
)
cal_spec <- deficiency_calibrate(delta, spec_equal)
stopifnot(abs(cal_spec$d - cal$d) < 1e-10)
stopifnot(abs(cal_spec$escape_mass - delta) < 1e-6)

eps_star <- (1 + 0.5)^(-q / 2)  # illustrative closure floor
const <- .c05_minorization_constants(eps_star, cal_spec$d, spec_equal)

stopifnot(abs(const$eps_d - eps_star * exp(-cal_spec$d)) < 1e-15)
stopifnot(abs(const$Q_mass_lb - .c05_weighted_chisq_cdf(2 * cal_spec$d, spec_equal$kappa)) < 1e-12)
stopifnot(const$eps == const$eps_d * const$Q_mass_lb)

n_sweeps <- ceiling(log(tol - delta) / log1p(-const$eps))
cat(sprintf(
  "q=%d delta=%.2f d=%.3f eps_star=%.4e eps=%.4e n=%d\n",
  q, delta, cal_spec$d, eps_star, const$eps, n_sweeps
))

cat("\n=== Checks B-D: Gaussian LMM on big_word_club (Ex_rLMM_reg path) ===\n")

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
mpl <- list(
  group.Sigma = ps$group.Sigma,
  group.dispersion = ps$group.dispersion,
  pop.prior_list = ps$pop.prior_list
)

mode_c05 <- population_mode(
  design = design,
  pfamily_list = pf,
  family = gaussian(),
  dispprior_list = dispprior_list,
  estep = "exact",
  icm_init = TRUE
)

mode_ref <- lmerb_posterior_mean(design, mpl)

max_gamma_gap <- max(vapply(
  design$groupef.names,
  function(k) max(abs(mode_c05$fixef[[k]] - mode_ref$fixef[[k]])),
  numeric(1L)
))
cat("max |gamma*_C05 - lmerb_posterior_mean|:", format(max_gamma_gap, digits = 3), "\n")
stopifnot(max_gamma_gap < 1e-8)

max_b_gap <- max(abs(mode_c05$b_mean - mode_ref$b_mean))
cat("max |b_mean gap|:", format(max_b_gap, digits = 3), "\n")
stopifnot(max_b_gap < 1e-8)

stopifnot(isTRUE(mode_c05$converged))
stopifnot(mode_c05$stationarity < 1e-8)
cat(
  "Gaussian EM iterations:", mode_c05$iterations,
  " ICM iterations:", mode_c05$icm$iterations, "\n"
)

cat("\n=== Check C: closure eps(gamma*) from tilde_J ===\n")

eps_obj <- epsilon_star(mode_c05, method = "closure")
manual <- exp(-0.5 * determinant(diag(mode_c05$q) + mode_c05$tilde_J,
                                 logarithm = TRUE)$modulus)
cat("eps_star:", signif(eps_obj$eps_star, 5),
    " manual:", signif(manual, 5), "\n")
stopifnot(abs(eps_obj$eps_star - manual) < 1e-10)
stopifnot(abs(eps_obj$eps_star - mode_c05$eps_star_closure) < 1e-15)

cat("\n=== Check C2: optimize matches closure on Gaussian LMM ===\n")

eps_opt <- epsilon_optimize(mode_c05, n = 2000L, mc_seed = 42L)
cat("optimize eps_star:", signif(eps_opt$eps_star, 5),
    " closure:", signif(eps_obj$eps_star, 5),
    " g_opt:", signif(eps_opt$g_opt, 5), "\n")
stopifnot(abs(eps_opt$eps_star - eps_obj$eps_star) < 1e-5)
stopifnot(isTRUE(eps_opt$certified))

cat("\n=== Check D: end-to-end certificate (Gaussian) ===\n")

cert <- certificate(
  design = design,
  pfamily_list = pf,
  family = gaussian(),
  delta = 0.01,
  dispprior_list = dispprior_list,
  tol = 0.05
)
print(cert)
stopifnot(cert$constants$eps > 0 && cert$constants$eps <= 1)
stopifnot(isTRUE(cert$certified))
stopifnot(identical(cert$epsilon$method, "closure"))
stopifnot(identical(cert$route, "spectrum_calibrated"))

cat("\nAll C05 epsilon scratch checks passed (Gaussian smoke; no binomial MC).\n")
cat("For full Gaussian parity see data-raw/_ex_c05_gaussian_parity.R\n")
