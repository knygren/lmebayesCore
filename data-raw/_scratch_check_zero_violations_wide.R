## TEMPORARY / SCRATCH -- follow-up to _scratch_check_zero_violations.R:
## same calibration, but with the package DEFAULT max_disp_perc_measurement
## (0.99, much wider truncation window) instead of Ex_13b's 0.8, to check
## whether the calibration search is a no-op ONLY because of the tight 0.8
## window, or unconditionally.

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
dat <- subset(dat, school_id %in% full_rank_schools)
dat$school_id <- droplevels(dat$school_id)

for (mdp in c(0.8, 0.9, 0.95, 0.99)) {
  ps <- Prior_Setup_GLMM(
    form_lmer,
    data            = dat,
    pwt             = 0.01,
    dispformula     = ~school_id,
    max_disp_perc_measurement = mdp,
    pwt_measurement = 0.1,
    alpha_target_measurement  = 0.01
  )
  calib <- ps$pwt_measurement_calibration
  cat(sprintf(
    "\nmax_disp_perc_measurement = %.2f:  pct_outside_before range [%.4g, %.4g], pct_outside_after range [%.4g, %.4g], groups w/ w_star>0: %d/%d, clipped: %d\n",
    mdp,
    min(calib$pct_outside_before), max(calib$pct_outside_before),
    min(calib$pct_outside_after), max(calib$pct_outside_after),
    sum(calib$w_star > 0), nrow(calib),
    sum(calib$clipped_at_ceiling)
  ))
}
