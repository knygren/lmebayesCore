## Scratch validation for refactored certificate pipeline (8 categories).
## Not part of the package test suite.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required to load lmebayesCore for this script.")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Suggested package 'bayesrules' is required for the Gaussian fixture.")
}

devtools::load_all("c:/Rpackages/lmebayesCore", quiet = TRUE)

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

cat("=== Category pipeline (explicit steps) ===\n")

mode <- population_mode(
  design = design, pfamily_list = pf, family = gaussian(),
  dispprior_list = dispprior_list, estep = "exact"
)
eps_obj <- epsilon_star(mode, method = "closure")
beta_mode <- beta_marginal_mode(
  design = design, pfamily_list = pf, family = gaussian(),
  dispprior_list = dispprior_list
)
beta_set <- beta_marginal_safe_set(delta_2 = 0.01, beta_mode = beta_mode)
floor_obj <- group_precision_floor(
  beta_mode = beta_mode, beta_set = beta_set, kappa_method = "none"
)
eigenvalues <- floor_coupling_eigenvalues(mode, floor_obj)
drift <- rosenthal_drift_constants(eigenvalues, display_mode = "sharp")

stopifnot(abs(max(eigenvalues$kappa_lb) - 0.9307937) < 1e-4)
stopifnot(abs(drift$lambda_lb - 0.9307937^2) < 1e-6)

opt <- optimal_rosenthal_tv_bound(
  mode = mode,
  eigenvalues = eigenvalues,
  eps = eps_obj$eps_star,
  k = 50L,
  delta_2 = 0.01,
  display_mode = "sharp"
)
stopifnot(opt$k == 50L)
stopifnot(opt$inner_bound > 0 && opt$inner_bound < 0.01)
stopifnot(opt$full_bound > opt$inner_bound)

cat("  k=50 inner:", signif(opt$inner_bound, 4),
    " full:", signif(opt$full_bound, 4),
    " kappa_max:", signif(eigenvalues$kappa_max_lb, 4), "\n")

cat("\n=== Orchestrator gamma_beta_tv_certificate ===\n")

cert <- gamma_beta_tv_certificate(
  design = design,
  pfamily_list = pf,
  family = gaussian(),
  delta_2 = 0.01,
  dispprior_list = dispprior_list,
  k = 50L,
  kappa_method = "none",
  estep = "mc",
  n = 2000L,
  mc_seed = 42L
)

stopifnot(inherits(cert, "gamma_beta_tv_certificate"))
stopifnot(inherits(cert$beta_mode, "beta_marginal_mode"))
stopifnot(is.null(cert$beta_set$Gamma_lb))
stopifnot(!is.null(cert$floor$Gamma_lb))
stopifnot(isTRUE(cert$certified$gamma_em))
stopifnot(isTRUE(cert$mode$restricted))
stopifnot(cert$rosenthal$bound > 0 && cert$rosenthal$bound < 1)

print(cert)

cat("\n=== Invert for inner tol = 0.01 ===\n")
cert_tol <- gamma_beta_tv_certificate(
  design = design, pfamily_list = pf, family = gaussian(),
  delta_2 = 0.01, dispprior_list = dispprior_list,
  inner_tol = 0.01, kappa_method = "none",
  estep = "mc", n = 2000L, mc_seed = 42L
)
cat("  k needed:", cert_tol$optimal$k,
    " inner:", signif(cert_tol$optimal$inner_bound, 4), "\n")

cat("\n=== Gaussian exact vs simulation (big_word_club) ===\n")

mc_n <- 10000L
mc_seed <- 42L

mode_exact <- population_mode(
  design = design, pfamily_list = pf, family = gaussian(),
  dispprior_list = dispprior_list, estep = "exact"
)
mode_mc <- population_mode(
  design = design, pfamily_list = pf, family = gaussian(),
  dispprior_list = dispprior_list, estep = "mc", n = mc_n, mc_seed = mc_seed
)

fixef_delta <- sqrt(sum(unlist(Map(function(a, b) a - b, mode_exact$fixef, mode_mc$fixef)))^2)
b_delta <- max(abs(mode_exact$b_mean - mode_mc$b_mean))
kappa_delta <- abs(mode_exact$kappa - mode_mc$kappa)
eps_closure_exact <- epsilon_star(mode_exact, method = "closure")
eps_closure_mc <- epsilon_star(mode_mc, method = "closure")
eps_opt_mc <- epsilon_star(mode_mc, method = "optimize", n = mc_n)

cat("  population_mode fixef L2 delta (exact vs MC n=", mc_n, "): ",
    signif(fixef_delta, 4), "\n", sep = "")
cat("  b_mean max |delta|:", signif(b_delta, 4), "\n")
cat("  kappa_max delta:", signif(kappa_delta, 6),
    " (rel ", signif(100 * kappa_delta / mode_exact$kappa, 4), "%)\n", sep = "")
cat("  eps_star closure (exact EM):", signif(eps_closure_exact$eps_star, 6), "\n")
cat("  eps_star closure (MC EM):   ", signif(eps_closure_mc$eps_star, 6),
    " rel bias vs exact: ",
    signif(100 * (eps_closure_mc$eps_star - eps_closure_exact$eps_star) /
             eps_closure_exact$eps_star, 4), "%\n", sep = "")
cat("  eps_star optimize (MC EM):  ", signif(eps_opt_mc$eps_star, 6),
    " rel vs closure-exact: ",
    signif(100 * (eps_opt_mc$eps_star - eps_closure_exact$eps_star) /
             eps_closure_exact$eps_star, 4), "%\n", sep = "")
cat("  eps optimize certified:", eps_opt_mc$certified,
    " mc_delta_floor:", signif(mode_mc$mc_delta_floor, 4), "\n")

estep_exact <- group_effects_conditional_mean(
  mode = mode_exact, estep = "exact"
)
estep_mc <- group_effects_conditional_mean(
  mode = mode_mc, estep = "mc", n = mc_n, mc_seed = mc_seed
)
estep_mc_at_exact <- group_effects_conditional_mean(
  mode = mode_exact, estep = "mc", n = mc_n, mc_seed = mc_seed
)
cat("  group_effects_conditional_mean b_mean max |delta|",
    " (exact mode+exact estep vs exact mode+MC): ",
    signif(max(abs(estep_exact$b_mean - estep_mc_at_exact$b_mean)), 4), "\n")

prior_list <- .c05_prior_list_from_fixef(
  design = design,
  fixef = mode_exact$fixef,
  p11 = mode_exact$p11,
  measurement_prior_list = mode_exact$measurement_prior_list,
  group_levels = mode_exact$p11$group_levels
)
draw_safe <- rNormal_reg_group_safe(
  n = 1L,
  y = design$y,
  x = design$D,
  group = design$group,
  prior_list = prior_list,
  beta_set = beta_set,
  mc_seed = mc_seed
)
draw_unsafe <- rNormal_reg_group(
  n = 1L,
  y = design$y,
  x = design$D,
  group = design$group,
  prior_list = prior_list
)
cat("  rNormal_reg_group_safe accept rate (1 draw, n=5000 target): ")
draw_safe_batch <- rNormal_reg_group_safe(
  n = 5000L,
  y = design$y,
  x = design$D,
  group = design$group,
  prior_list = prior_list,
  beta_set = beta_set,
  mc_seed = mc_seed
)
cat(signif(draw_safe_batch$accept_rate, 4),
    " (", draw_safe_batch$n_tried, " tries)\n", sep = "")
cat("  safe vs exact b_mean max |delta| (5000 accepted): ",
    signif(max(abs(draw_safe_batch$b_mean - estep_exact$b_mean)), 4), "\n", sep = "")

mode_restricted <- population_mode(
  design = design, pfamily_list = pf, family = gaussian(),
  dispprior_list = dispprior_list, estep = "exact",
  beta_set = beta_set, n = 5000L, mc_seed = mc_seed
)
cat("  population_mode with beta_set (5000 accepted draws): accept_rate=",
    signif(mode_restricted$accept_rate, 4),
    " fixef L2 delta vs exact:", signif(
      sqrt(sum(unlist(Map(function(a, b) a - b,
                          mode_exact$fixef, mode_restricted$fixef)))^2),
      4), "\n", sep = "")

cat("\nAll refactored certificate checks passed.\n")
