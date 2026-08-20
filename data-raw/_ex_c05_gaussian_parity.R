# Gaussian-only parity: population_mode vs lmerb / rlmerb / rglmerb.
# Complex Gaussian LMM on bayesrules::big_word_club (inst/examples/Ex_rLMM_reg.R):
#   random intercept + slopes, population mean slopes, and cross-level
#   moderation of a random slope (free_reduced_lunch:distracted_a1 in W).
# Not part of the package test suite.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools required to load lmebayesCore for this script.")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop(
    "Suggested package 'bayesrules' is required ",
    "(same data path as inst/examples/Ex_rLMM_reg.R)."
  )
}
devtools::load_all(".", quiet = TRUE)

.max_fixef_gap <- function(a, b, re_names) {
  max(vapply(
    re_names,
    function(k) max(abs(a[[k]] - b[[k]])),
    numeric(1L)
  ))
}

.max_b_gap <- function(a, b) {
  b <- b[rownames(a), colnames(a), drop = FALSE]
  max(abs(a - b))
}

data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
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
re_names <- design$groupef.names

cat("=== Design: random slopes + population predictors of slopes ===\n")
for (k in re_names) {
  cat("W[[", k, "]] cols:", paste(colnames(design$W[[k]]), collapse = ", "), "\n")
}
stopifnot(
  "(Intercept)" %in% re_names,
  "distracted_ppvt" %in% re_names,
  "distracted_a1" %in% re_names,
  ncol(design$W[["distracted_a1"]]) >= 2L,
  "free_reduced_lunch" %in% colnames(design$W[["distracted_a1"]])
)

mpl_ps <- list(
  group.Sigma = ps$group.Sigma,
  group.dispersion = ps$group.dispersion,
  pop.prior_list = ps$pop.prior_list
)

# rlmerb / population_mode resolve measurement priors from pfamily_list
# (Block~1 tau^2 plug-ins), not from Prior_Setup_GLMM$group.Sigma directly.
.priors_from_pfamily_list <- getFromNamespace(
  ".priors_from_pfamily_list", "lmebayesCore"
)
prior_pack <- .priors_from_pfamily_list(
  pfamily_list = pf,
  group.dispersion = ps$group.dispersion,
  design = design,
  family = gaussian(),
  fn_name = "parity"
)
mpl <- list(
  group.Sigma = prior_pack$group.Sigma,
  group.dispersion = prior_pack$group.dispersion,
  pop.prior_list = prior_pack$pop.prior_list
)

cat("\n=== Reference: lmerb_posterior_mean (pfamily_list measurement prior) ===\n")
ref <- lmerb_posterior_mean(design, mpl)
print(ref$fixef)

cat(
  "\n(note: Prior_Setup group.Sigma vs pfamily path max |diff|:",
  format(max(abs(mpl_ps$group.Sigma - mpl$group.Sigma)), digits = 3), ")\n"
)

cat("\n=== population_mode (C05 EM, estep = exact) ===\n")
mode_c05 <- population_mode(
  design = design,
  pfamily_list = pf,
  family = gaussian(),
  dispprior_list = dispprior_list,
  estep = "exact",
  icm_init = TRUE
)
gap_gamma <- .max_fixef_gap(mode_c05$fixef, ref$fixef, re_names)
gap_b <- .max_b_gap(mode_c05$b_mean, ref$b_mean)
cat(
  "max |gamma - ref|:", format(gap_gamma, digits = 3),
  "  max |b - ref|:", format(gap_b, digits = 3),
  "  EM iters:", mode_c05$iterations,
  "  q =", mode_c05$q, "\n"
)
stopifnot(gap_gamma < 1e-8, gap_b < 1e-8)
stopifnot(isTRUE(mode_c05$converged))

cat("\n=== rlmerb(..., simulate = FALSE) popef.mode (= exact mean) ===\n")
pt_lmm <- rlmerb(
  design = design,
  pfamily_list = pf,
  dispprior_list = dispprior_list,
  simulate = FALSE,
  verbose = FALSE,
  print_icm_table = FALSE
)
ref_rlmerb <- lmerb_posterior_mean(design, list(
  group.Sigma = pt_lmm$prior$group.Sigma,
  group.dispersion = pt_lmm$prior$group.dispersion,
  pop.prior_list = pt_lmm$prior$pop.prior_list
))
gap_internal <- .max_fixef_gap(ref_rlmerb$fixef, pt_lmm$popef.mode, re_names)
cat("rlmerb internal consistency |lmerb - popef.mode|:", format(gap_internal, digits = 3), "\n")
stopifnot(gap_internal < 1e-10)

gap_rlmerb <- .max_fixef_gap(pt_lmm$popef.mode, mode_c05$fixef, re_names)
gap_rlmerb_b <- .max_b_gap(pt_lmm$groupef.mode, mode_c05$b_mean)
cat(
  "max |popef.mode - population_mode|:", format(gap_rlmerb, digits = 3),
  "  max |groupef.mode - b_mean|:", format(gap_rlmerb_b, digits = 3), "\n"
)
stopifnot(gap_rlmerb < 1e-8, gap_rlmerb_b < 1e-8)

cat("\n=== rglmerb(..., simulate = FALSE, family = gaussian()) ===\n")
pt_glm <- rglmerb(
  design = design,
  pfamily_list = pf,
  family = gaussian(),
  dispprior_list = dispprior_list,
  simulate = FALSE,
  verbose = FALSE
)
gap_rglmerb <- .max_fixef_gap(pt_glm$popef.mode, mode_c05$fixef, re_names)
gap_rglmerb_b <- .max_b_gap(pt_glm$groupef.mode, mode_c05$b_mean)
cat(
  "max |popef.mode - population_mode|:", format(gap_rglmerb, digits = 3),
  "  max |groupef.mode - b_mean|:", format(gap_rglmerb_b, digits = 3), "\n"
)
stopifnot(gap_rglmerb < 1e-8, gap_rglmerb_b < 1e-8)

cat("\n=== rlmerb iid draws: popef.means vs exact mean ===\n")
set.seed(20260818)
fit <- rlmerb(
  n = 200L,
  design = design,
  pfamily_list = pf,
  dispprior_list = dispprior_list,
  simulate = TRUE,
  verbose = FALSE,
  progbar = FALSE,
  print_icm_table = FALSE
)
stopifnot(identical(fit$convergence_info$draw_engine, "rLMMNormal_joint_iid"))
gap_means <- .max_fixef_gap(fit$popef.means, mode_c05$fixef, re_names)
cat(
  "max |popef.means - population_mode| (n=200 iid):",
  format(gap_means, digits = 3), "\n"
)
stopifnot(gap_means < 0.35)

cat("\n=== C05 closure / optimize on big_word_club ===\n")
eps_closure <- epsilon_star(mode_c05, method = "closure")
eps_opt <- epsilon_optimize(mode_c05)
cat(
  "eps_star closure:", signif(eps_closure$eps_star, 5),
  " optimize:", signif(eps_opt$eps_star, 5), "\n"
)
stopifnot(abs(eps_opt$eps_star - eps_closure$eps_star) < 1e-5)
stopifnot(isTRUE(eps_opt$certified))

cat("\nAll Gaussian parity checks passed.\n")
