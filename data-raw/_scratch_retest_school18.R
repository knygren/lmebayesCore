## Retest: does school_id "18" still trigger the "Sign violation: UB2 < 0"
## Block~1 ING envelope failure under the CURRENT prior specification
## (Part VI default + n_j/2-corrected disp bounds + alpha_target_measurement
## calibration) and the CURRENT glmbayesCore (root-finding UB2_Min_j fix,
## data-raw/README_ub2_rootfinding_fix.md)? Mirrors Ex_13b's Section 1/2/5
## but keeps school 18 in (still drops "2" to isolate 18 specifically).
devtools::load_all(".", quiet = TRUE)

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

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

design_all <- model_setup(form_lmer, data = dat)
full_rank_schools <- names(design_all$re_rank)[design_all$re_rank]
cat(sprintf(
  "Full-rank filter: %d of %d schools kept\n",
  length(full_rank_schools), length(design_all$re_rank)
))
dat <- subset(dat, school_id %in% full_rank_schools)
dat$school_id <- droplevels(dat$school_id)

## Keep BOTH "18" and "2" in this time (no TEMP drop at all).
cat("school_id '18' present:", "18" %in% levels(dat$school_id), "\n")
cat("school_id '2' present:", "2" %in% levels(dat$school_id), "\n")
cat("n_j for school 18:", sum(dat$school_id == "18"), "\n")
cat("n_j for school 2:", sum(dat$school_id == "2"), "\n")

design <- model_setup(form_lmer, data = dat)
stopifnot(all(design$re_rank))

ps <- Prior_Setup_GLMM(
  form_lmer,
  data            = dat,
  pwt             = 0.01,
  dispformula     = ~school_id,
  max_disp_perc_measurement = 0.8,
  pwt_measurement = 0.1,
  alpha_target_measurement  = 0.01
)

grp <- design$group
attr(grp, "group_name") <- design$group_name
group_levels <- levels(grp)
re_names <- design$re_coef_names
p_re <- length(re_names)
pf <- pfamily_list(ps)

disp_pf_list <- dGamma_list(ps, max_disp_perc_measurement = 0.8)
shape_group      <- stats::setNames(numeric(length(group_levels)), group_levels)
rate_group       <- stats::setNames(numeric(length(group_levels)), group_levels)
disp_lower_group <- stats::setNames(numeric(length(group_levels)), group_levels)
disp_upper_group <- stats::setNames(numeric(length(group_levels)), group_levels)
for (lev in group_levels) {
  pl <- disp_pf_list[[lev]]$prior_list
  shape_group[[lev]]      <- pl$shape[1L]
  rate_group[[lev]]       <- pl$rate[1L]
  disp_lower_group[[lev]] <- pl$disp_lower
  disp_upper_group[[lev]] <- pl$disp_upper
}
cat(sprintf(
  "school 18: sigma2_hat = %.4f, disp_lower = %.4f, disp_upper = %.4f\n",
  rate_group[["18"]] / (shape_group[["18"]] - 1),
  disp_lower_group[["18"]], disp_upper_group[["18"]]
))

prior_list <- list(
  mu               = matrix(0, nrow = p_re, ncol = 1L, dimnames = list(re_names, NULL)),
  Sigma            = as.matrix(ps$Sigma_ranef),
  shape_group      = shape_group,
  rate_group       = rate_group,
  disp_lower_group = disp_lower_group,
  disp_upper_group = disp_upper_group
)

cat("\n--- Running sampler (n = 200, verbose = FALSE, all groups incl. 18 AND 2) ---\n\n")
fit <- rLMMindepNormalGamma_reg_known_vcov(
  n            = 200L,
  y            = design$y,
  D            = design$D,
  group        = grp,
  W            = design$W,
  pfamily_list = pf,
  dispprior_list = prior_list,
  progbar      = FALSE,
  verbose      = FALSE
)

cat("\n=== Sampler completed with NO error ===\n")
cat("m_convergence =", fit$m_convergence, "\n")
cat("all dispersion_ranef finite/positive:",
    all(is.finite(fit$group.dispersion)) && all(fit$group.dispersion > 0), "\n")
cat("school 18 dispersion_ranef.mean:", fit$group.dispersion.mean[["18"]], "\n")
cat("school 2 dispersion_ranef.mean:", fit$group.dispersion.mean[["2"]], "\n")
