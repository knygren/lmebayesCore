devtools::load_all(".", quiet = TRUE)
data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) & !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)
form <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)
design0 <- model_setup(form, data = dat)
dat <- subset(dat, school_id %in% names(design0$groupef.rank)[design0$groupef.rank])
dat$school_id <- droplevels(dat$school_id)
design <- model_setup(form, data = dat)
ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)
pf <- pfamily_list(ps)
dispprior_list <- list(dispersion = ps$group.dispersion)

.priors_from_pfamily_list <- getFromNamespace(".priors_from_pfamily_list", "lmebayesCore")
prior_pack <- .priors_from_pfamily_list(pf, ps$group.dispersion, design, gaussian(), "cmp")
mpl <- list(
  group.Sigma = prior_pack$group.Sigma,
  group.dispersion = prior_pack$group.dispersion,
  pop.prior_list = prior_pack$pop.prior_list
)

mode <- population_mode(
  design, pf, gaussian(), dispprior_list,
  estep = "exact", icm_init = TRUE
)

block1 <- list(
  Sigma = prior_pack$group.Sigma,
  dispersion = prior_pack$group.dispersion
)
rate <- two_block_rate_from_pfamily_list(
  x = design$D,
  group = design$group,
  x_hyper = design$W,
  prior_list_block1 = block1,
  pfamily_list = pf,
  family = gaussian()
)

cat("=== C05 population_mode ===\n")
cat("kappa (max coupling spectrum):", signif(mode$kappa, 5), "\n")
cat("kappa_spectrum:", paste(signif(mode$kappa_spectrum, 4), collapse = ", "), "\n")
cat("eps_star_closure:", signif(mode$eps_star_closure, 6), "\n")

cat("\n=== two_block_rate (Remark 8) ===\n")
cat("lambda*:", signif(rate$lambda_star, 5), "\n")
cat("1 - lambda*:", signif(1 - rate$lambda_star, 5), "\n")
cat("top eigenvalues:", paste(signif(rate$eigenvalues, 4), collapse = ", "), "\n")

q <- mode$q
lam <- rate$lambda_star
rough <- (1 + lam)^(-q / 2)
cat("\nRough (1+lambda*)^(-q/2):", signif(rough, 5), "\n")

const <- epsilon(mode$eps_star_closure, delta = 0.01, mode = mode)
cat("\n=== Full certificate pipeline (delta=0.01) ===\n")
cat("eps_star (Stage 2):", signif(const$eps_star, 6), "\n")
cat("d (spectrum):", signif(const$d, 4), "  exp(-d):", signif(exp(-const$d), 4), "\n")
cat("eps_d = exp(-d)*eps_star:", signif(const$eps_d, 4), "\n")
cat("Q_mass_lb:", signif(const$Q_mass_lb, 4), "\n")
cat("eps (Theorem 2 multiplier):", signif(const$eps, 6), "\n")
