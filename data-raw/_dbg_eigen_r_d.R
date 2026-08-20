# Report coupling eigenvalues and r/d calibration (B-spectrum route).
devtools::load_all(".", quiet = TRUE)
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
form <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)
design0 <- model_setup(form, data = dat_g)
dat_g <- subset(dat_g, school_id %in% names(design0$groupef.rank)[design0$groupef.rank])
dat_g$school_id <- droplevels(dat_g$school_id)
design <- model_setup(form, data = dat_g)
ps <- Prior_Setup_GLMM(form, data = dat_g, pop.pwt = 0.01)
pf <- pfamily_list(ps)
mode <- population_mode(
  design, pf, gaussian(),
  dispprior_list = list(dispersion = ps$group.dispersion),
  estep = "exact", icm_init = TRUE
)
spec <- deficiency_spectrum(mode)
delta <- 0.01
cal <- deficiency_calibrate(delta, spec)
eps_star <- epsilon_star(mode, method = "closure")$eps_star
const <- epsilon(eps_star, delta, mode)

cat("=== A eigenvalues kappa_i = eig(P11^{-1/2} S P11^{-1/2}) ===\n")
print(signif(spec$kappa, 5))
cat("\n=== B weights w_i = kappa_i/(1-kappa_i) = eig(B) ===\n")
print(signif(spec$weights, 5))
cat("\nsum(kappa) =", signif(sum(spec$kappa), 4),
    "  kappa_max =", signif(spec$kappa_max, 5),
    "  rho =", signif(spec$rho, 4), "\n")
cat("\n=== Calibration: Pr(sum w_i Z_i^2 > r^2) = delta, d = r^2/2 ===\n")
cat("delta =", delta, "\n")
cat("r^2 =", signif(cal$r2, 5), "  r =", signif(cal$r, 5),
    "  d =", signif(cal$d, 5), "\n")
cat("escape_mass =", signif(cal$escape_mass, 4), "\n")
cat("Q_mass_lb (kappa, r^2) =",
    signif(.c05_weighted_chisq_cdf(cal$r2, spec$kappa), 4), "\n")
cat("eps(gamma*) =", signif(eps_star, 5),
    "  eps_d =", signif(const$eps_d, 4),
    "  eps =", signif(const$eps, 4), "\n")

cat("\n=== Limit checks (synthetic q=1; eps_star = (1+kappa)^{-1/2} from closure) ===\n")
cat("Note: eps(gamma*) falls as kappa rises (stronger coupling / curvature).\n\n")
cat(sprintf(
  "%-7s %-9s %-7s %-7s %-11s %-11s %-11s\n",
  "kappa", "w", "r", "d", "eps_star", "eps_d", "eps"
))
for (k in c(0.05, 0.5, 0.95)) {
  w1 <- k / (1 - k)
  r1 <- .c05_r_from_delta(delta, w1)
  eps_star_k <- (1 + k)^(-1 / 2)
  Q_lb <- stats::pchisq(r1$r2 / k, df = 1)
  eps_d <- exp(-r1$d) * eps_star_k
  eps <- eps_d * Q_lb
  cat(sprintf(
    "%-7.2f %-9.2f %-7.3f %-7.3f %-11.5f %-11.4e %-11.4e\n",
    k, w1, r1$r, r1$d, eps_star_k, eps_d, eps
  ))
}

cat("\n=== Full spectrum (big_word_club) ===\n")
cat(sprintf(
  "%-7s %-9s %-7s %-7s %-11s %-11s %-11s\n",
  "kappa", "w", "r", "d", "eps_star", "eps_d", "eps"
))
cat(sprintf(
  "multi   multi     %-7.3f %-7.3f %-11.5f %-11.4e %-11.4e\n",
  const$r, const$d, eps_star, const$eps_d, const$eps
))
