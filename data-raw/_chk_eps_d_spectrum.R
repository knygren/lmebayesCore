# Scratch checks for spectrum-based deficiency calibration (C05 §4A).
# Not part of the package test suite.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required to load lmebayesCore for this script.")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Suggested package 'bayesrules' is required.")
}
devtools::load_all(".", quiet = TRUE)

cat("=== Check 1: equal-kappa B-spectrum ===\n")
q <- 7L
delta <- 0.01
k <- 0.5
spec_equal <- list(
  kappa = rep(k, q),
  weights = rep(k / (1 - k), q),
  kappa_max = k,
  q = q,
  method = "symmetrized_coupling"
)
cal <- deficiency_calibrate(delta, spec_equal)
r2_ref <- (k / (1 - k)) * stats::qchisq(1 - delta, df = q)
stopifnot(abs(cal$r2 - r2_ref) < 1e-8)
stopifnot(abs(cal$d - 0.5 * cal$r2) < 1e-15)
stopifnot(abs(cal$escape_mass - delta) < 1e-6)
cat("  r =", signif(cal$r, 5), " d =", signif(cal$d, 5), "\n")

cat("\n=== Check 2: big_word_club Gaussian closure ===\n")
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

mode <- population_mode(
  design = design,
  pfamily_list = pf,
  family = gaussian(),
  dispprior_list = dispprior_list,
  estep = "exact",
  icm_init = TRUE
)
stopifnot(isTRUE(mode$converged))

spec <- deficiency_spectrum(mode)
eps_star <- epsilon_star(mode, method = "closure")$eps_star
const <- epsilon(eps_star, delta, mode)

cat("  kappa (A):", paste(signif(spec$kappa, 4), collapse = ", "), "\n")
cat("  weights w (B):", paste(signif(spec$weights, 4), collapse = ", "), "\n")
cat("  r =", signif(const$r, 5), " d =", signif(const$d, 5), "\n")
cat("  escape_mass:", signif(const$cal$escape_mass, 4), "\n")

stopifnot(abs(const$cal$escape_mass - delta) < 1e-4)
stopifnot(abs(const$d - 0.5 * const$r2) < 1e-12)
stopifnot(abs(const$eps_d - eps_star * exp(-const$d)) < 1e-15)
stopifnot(const$eps == const$eps_d * const$Q_mass_lb)
stopifnot(identical(const$route, "spectrum_calibrated"))
stopifnot(all(abs(spec$weights - spec$kappa / (1 - spec$kappa)) < 1e-12))

cat("\nAll spectrum deficiency calibration checks passed.\n")
