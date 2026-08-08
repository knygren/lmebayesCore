## Scratch: quick end-to-end smoke test of the new split-support revised TV
## bound (inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 17,
## data-raw/_scratch_tv_bound_revised.R) against a small
## rLMMindepNormalGamma_reg_known_vcov() fit -- small n, verbose = FALSE, no
## group drop (Ex_13b-like fixture). Speed, not final numbers, is the point.
devtools::load_all(".", quiet = TRUE)

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
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
form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

design_all <- model_setup(form_lmer, data = dat)
full_rank_schools <- names(design_all$re_rank)[design_all$re_rank]
dat <- subset(dat, school_id %in% full_rank_schools)
dat$school_id <- droplevels(dat$school_id)
dat <- subset(dat, !as.character(school_id) %in% c("18", "2"))
dat$school_id <- droplevels(dat$school_id)

design <- model_setup(form_lmer, data = dat)
stopifnot(all(design$re_rank))

ps <- Prior_Setup_GLMM(
  form_lmer, data = dat, pwt = 0.01, dispformula = ~school_id,
  max_disp_perc_measurement = 0.8, pwt_measurement = 0.1
)
pf <- pfamily_list(ps)
disp_pf_list <- dGamma_list(ps, max_disp_perc_measurement = 0.8)

grp <- design$group
attr(grp, "group_name") <- design$group_name
group_levels <- levels(grp)
re_names     <- design$re_coef_names
p_re         <- length(re_names)

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

prior_list <- list(
  mu               = matrix(0, nrow = p_re, ncol = 1L, dimnames = list(re_names, NULL)),
  Sigma            = as.matrix(ps$Sigma_ranef),
  shape_group      = shape_group,
  rate_group       = rate_group,
  disp_lower_group = disp_lower_group,
  disp_upper_group = disp_upper_group
)

cat(sprintf(
  "=== Running rLMMindepNormalGamma_reg_known_vcov(n = 60) on %d groups ===\n",
  length(group_levels)
))

fit <- rLMMindepNormalGamma_reg_known_vcov(
  n            = 60L,
  y            = design$y,
  D            = design$D,
  group        = grp,
  W            = design$W,
  pfamily_list = pf,
  dispprior_list = prior_list,
  progbar      = FALSE,
  verbose      = FALSE
)
n_draws <- nrow(fit$popef[[re_names[1L]]])
stopifnot(n_draws == 60L)

cat(sprintf(
  "m_convergence = %s, lambda_star_marginal = %s, marginal_rate_valid = %s, tv_tol = %s\n",
  format(fit$convergence_info$m_convergence),
  format(fit$convergence_info$lambda_star_marginal),
  format(fit$convergence_info$marginal_rate_valid),
  format(fit$convergence_info$tv_tol)
))

## Same Section 7d recipe as Ex_13b/Ex_13c.
source("data-raw/_scratch_lambda_star_marginal_over_draws.R")
prior_list_block1_rate <- list(
  Sigma      = as.matrix(ps$Sigma_ranef),
  dispersion = disp_upper_group[group_levels]
)
prior_list_block2_rate <- lapply(pf, function(pfk) {
  pl <- pfk$prior_list
  list(mu = pl$mu, Sigma = pl$Sigma, dispersion = pl$dispersion)
})
inp_marg <- lmebayesCore:::.two_block_rate_inputs(
  x = design$D, block = grp, x_hyper = design$W,
  prior_list_block1 = prior_list_block1_rate,
  prior_list_block2 = prior_list_block2_rate
)
blocks_marg <- lmebayesCore:::.two_block_S_P11(inp_marg)
group_setup_marg <- .tmp_marginal_group_setup(
  D = design$D, y = design$y, group = grp, group_levels = group_levels,
  re_coef_names = re_names, shape_group = shape_group, rate_group = rate_group
)
res_marg <- .tmp_lambda_star_marginal_over_draws(
  fit = fit, n_draws = n_draws, y = design$y,
  group_name = design$group_name, group_setup = group_setup_marg,
  inp = inp_marg, blocks = blocks_marg
)
.tmp_print_marginal_over_draws_summary(res_marg)

## New Section 7e: the revised TV bound under test.
source("data-raw/_scratch_tv_bound_revised.R")
lambda_star_used <- if (isTRUE(fit$convergence_info$marginal_rate_valid)) {
  fit$convergence_info$lambda_star_marginal
} else {
  fit$convergence_info$lambda_star_upper
}
tv_revised <- .tmp_tv_bound_revised(
  res_marg, lambda_star_used, tv_tol = fit$convergence_info$tv_tol
)
.tmp_print_tv_bound_revised(tv_revised)

stopifnot(is.data.frame(tv_revised), nrow(tv_revised) == 1L)
stopifnot(all(is.finite(unlist(tv_revised))))
stopifnot(tv_revised$p_hat >= 0, tv_revised$p_hat <= 1)
stopifnot(abs(tv_revised$tv_revised - (tv_revised$tv_tol + 2 * tv_revised$p_hat)) < 1e-12)
stopifnot(abs(tv_revised$tv_revised_loose - (tv_revised$tv_tol + tv_revised$p_hat)) < 1e-12)
stopifnot(tv_revised$tv_revised >= tv_revised$tv_revised_loose)
stopifnot(tv_revised$tv_revised >= tv_revised$tv_tol)

cat("\nOK: .tmp_tv_bound_revised() runs end-to-end and is internally consistent.\n")
