## Ex_21: Prior_Setup_lmebayes only (no lmerb sampling)
pkgload::load_all(".", quiet = TRUE)

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("This script requires the 'lme4' package.", call. = FALSE)
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

design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup ===\n\n")
print(design)

ps <- Prior_Setup_lmebayes(
  form_lmer,
  data           = dat,
  pwt            = 0.01,
  pwt_dispersion = 0.2
)

cat("\n=== print(Prior_Setup_lmebayes) ===\n\n")
print(ps)

cat("\n=== DONE: Ex_21 prior setup only ===\n")
