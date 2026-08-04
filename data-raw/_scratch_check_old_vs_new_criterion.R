## TEMPORARY / SCRATCH -- reconciles the "several groups above threshold"
## result from the earlier (UNTRUNCATED, Section 16.3) pwt_measurement
## search against the current (EXACT/TRUNCATED, Section 16.6) one, on the
## SAME ps/design, to show both are internally consistent -- they are
## answering two different questions.

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

design <- model_setup(form_lmer, data = dat)
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
re_names <- design$re_coef_names

source("data-raw/_scratch_group_pwt_measurement_noncentral.R")
tab_pwt <- .tmp_group_pwt_measurement_noncentral(
  ps = ps, design = design, group = grp, group_levels = group_levels,
  re_coef_names = re_names
)
cat("\n=== OLD untruncated (Section 16.3) vs NEW exact/truncated (Section 16.6) criterion ===\n\n")
print(
  tab_pwt[order(-tab_pwt$w_star_noncentral_untrunc),
          c("group", "n_j", "pct_outside_at_current", "pct_outside_at_current_exact",
            "w_star_noncentral_untrunc", "w_star_exact")],
  row.names = FALSE, digits = 4
)
