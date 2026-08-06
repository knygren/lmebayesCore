## TEMPORARY / SCRATCH -- follow-up: confirm the calibration search itself is
## not broken (i.e. it *can* detect violations) by re-running it with a
## deliberately inadequate floor (pwt_measurement = 1e-4, near w_lo) instead
## of the real 0.1 default, and checking whether pct_outside_before becomes
## nonzero and uniroot() actually engages for at least some groups.

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

ps <- Prior_Setup_GLMM(
  form_lmer,
  data            = dat,
  pwt             = 0.01,
  dispformula     = ~school_id,
  max_disp_perc_measurement = 0.8,
  pwt_measurement = 1e-4,   ## deliberately tiny floor
  alpha_target_measurement  = 0.01
)
calib <- ps$pwt_measurement_calibration
cat("\n=== calibration with pwt_measurement floor = 1e-4 ===\n\n")
print(calib, row.names = FALSE, digits = 6)
cat(sprintf(
  "\nGroups w/ w_star > 0 (search actually engaged): %d / %d\n",
  sum(calib$w_star > 0), nrow(calib)
))
cat(sprintf(
  "Range pct_outside_before: [%.4g, %.4g]\n",
  min(calib$pct_outside_before), max(calib$pct_outside_before)
))
