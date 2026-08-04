## Scratch smoke test: verify the new .two_block_combine_rate_envelopes()
## path (rank-matched max of the disp_upper plug-in envelope and the
## Omega-marginalized envelope) fires and produces sane output. Mirrors
## demo/Ex_13b_rLMM_estimated_dispersion_known_vcov_BigWordClub_PartVI.R
## Sections 1-5 but with a tiny n so it runs in seconds.
##
## Not a permanent test -- throwaway, per .cursor/rules/testthat-and-scratch-
## code.mdc.

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

design <- model_setup(form_lmer, data = dat)
stopifnot(all(design$re_rank))

ps <- Prior_Setup_lmebayes(
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

pf <- pfamily_list(ps)
disp_pf_list <- dGamma_list(ps, max_disp_perc_measurement = 0.8)

set.seed(1)
fit <- rLMMindepNormalGamma_reg_known_vcov(
  n            = 20L,
  y            = design$y,
  D            = design$D,
  group        = grp,
  W            = design$W,
  prior_list   = disp_pf_list,
  pfamily_list = pf,
  progbar      = FALSE,
  verbose      = TRUE
)

ci <- fit$convergence_info
cat("\n=== convergence_info summary ===\n\n")
cat(sprintf("marginal_rate_valid   = %s\n", ci$marginal_rate_valid))
cat(sprintf("lambda_star_upper     = %.4f\n", ci$lambda_star_upper))
cat(sprintf("lambda_star_marginal  = %.4f\n", ci$lambda_star_marginal))
cat(sprintf("lambda_star_combined  = %s\n",
            if (is.null(ci$lambda_star_combined)) "NULL" else sprintf("%.4f", ci$lambda_star_combined)))
cat(sprintf("m_min_upper           = %s\n", ci$m_min_upper))
cat(sprintf("m_min_marginal        = %s\n", ci$m_min_marginal))
cat(sprintf("m_min_combined        = %s\n", if (is.null(ci$m_min_combined)) "NULL" else ci$m_min_combined))
cat(sprintf("m_convergence (used)  = %d\n", ci$m_convergence))

stopifnot(
  isTRUE(ci$marginal_rate_valid),
  !is.null(ci$lambda_star_combined),
  !is.null(ci$m_min_combined),
  ci$lambda_star_combined >= ci$lambda_star_upper - 1e-8,
  ci$lambda_star_combined >= ci$lambda_star_marginal - 1e-8,
  ci$m_convergence >= ci$m_min_upper,
  ci$m_convergence >= ci$m_min_marginal
)
cat("\nAll checks passed.\n")
