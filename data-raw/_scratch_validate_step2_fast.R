## TEMPORARY / SCRATCH -- fast (no sampler) validation of Step 2's exact-
## criterion pwt_measurement derivation. Not package code.

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

## Match Ex_13b's extra manual drop for Block~1 ING envelope failure.
temp_drop_schools <- c("18", "2")
drop <- intersect(temp_drop_schools, levels(dat$school_id))
if (length(drop)) {
  dat <- subset(dat, !as.character(school_id) %in% drop)
  dat$school_id <- droplevels(dat$school_id)
}

design <- model_setup(form_lmer, data = dat)
stopifnot(all(design$re_rank))

ps_flat <- Prior_Setup_GLMM(
  form_lmer,
  data            = dat,
  pwt             = 0.01,
  dispformula     = ~school_id,
  max_disp_perc_measurement = 0.8,
  pwt_measurement = 0.1
)

grp <- design$group
attr(grp, "group_name") <- design$group_name
group_levels <- levels(grp)
re_names     <- design$re_coef_names

source("data-raw/_scratch_group_pwt_measurement_noncentral.R")
tab_pwt <- .tmp_group_pwt_measurement_noncentral(
  ps = ps_flat, design = design, group = grp, group_levels = group_levels,
  re_coef_names = re_names
)
.tmp_print_group_pwt_measurement_noncentral(tab_pwt)

cat("\n\n=== Cross-check: does w_star_exact_floored actually bring pct_outside_exact below target? ===\n\n")

## Re-derive, for each group, pct_outside_exact AT the chosen w_star_exact
## (not w_star_exact_floored, since flooring can only help further) to
## confirm the uniroot() actually solved the intended equation.
verify <- do.call(rbind, lapply(seq_len(nrow(tab_pwt)), function(i) {
  lev <- tab_pwt$group[i]
  data.frame(
    group = lev,
    w_exact = tab_pwt$w_star_exact[i],
    clipped = tab_pwt$clipped_at_0.5_exact[i],
    pct_outside_now_untrunc = tab_pwt$pct_outside_at_current[i],
    pct_outside_now_exact = tab_pwt$pct_outside_at_current_exact[i]
  )
}))
print(verify, row.names = FALSE, digits = 4)

cat(sprintf(
  "\nMean pct_outside_at_current (untruncated) vs exact, across groups: %.3f%% vs %.3f%%\n",
  mean(tab_pwt$pct_outside_at_current), mean(tab_pwt$pct_outside_at_current_exact)
))
cat(sprintf(
  "Groups with exact pct_outside_at_current > 1%%: %d / %d (untruncated: %d / %d)\n",
  sum(tab_pwt$pct_outside_at_current_exact > 1), nrow(tab_pwt),
  sum(tab_pwt$pct_outside_at_current > 1), nrow(tab_pwt)
))

cat("\nDone.\n")
