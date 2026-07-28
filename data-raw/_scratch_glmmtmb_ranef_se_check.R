## Scratch: inspect the structure of glmmTMB::ranef(fit, condVar = TRUE) to
## find the right way to get per-group/per-coefficient conditional SEs for
## standardizing (z-scoring) the random effects.

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

fit <- glmmTMB::glmmTMB(form_lmer, data = dat, dispformula = ~school_id)

re <- glmmTMB::ranef(fit, condVar = TRUE)
str(re, max.level = 2)
cat("\nclass(re$cond$school_id):", class(re$cond$school_id), "\n")
cv <- attr(re$cond$school_id, "condVar")
cat("\nclass(condVar attr):", class(cv), "\n")
cat("dim(condVar attr):", paste(dim(cv), collapse = " x "), "\n")
print(cv[, , 1])

## Standardize: z_jk = ranef_jk / sqrt(condVar_kk,j)
re_df <- re$cond$school_id
se <- t(apply(cv, 3, function(m) sqrt(diag(m))))
colnames(se) <- colnames(re_df)
rownames(se) <- rownames(re_df)
z <- re_df / se

cat("\n=== Standardized (z-scored) random effects ===\n\n")
print(round(z, 3))

## Bonferroni-corrected threshold across all groups x RE coefficients.
n_tests <- length(z)
z_crit <- qnorm(1 - 0.05 / (2 * n_tests))
cat(sprintf(
  "\nBonferroni z threshold (n_tests = %d, alpha = 0.05): %.3f\n",
  n_tests, z_crit
))
flagged <- which(abs(as.matrix(z)) > z_crit)
cat("Flagged (group, coef) cells (|z| > threshold):\n")
if (length(flagged)) {
  idx <- arrayInd(flagged, dim(z))
  for (i in seq_len(nrow(idx))) {
    cat(sprintf("  %-4s  %-18s  z = %.3f\n",
                rownames(z)[idx[i, 1]], colnames(z)[idx[i, 2]],
                z[idx[i, 1], idx[i, 2]]))
  }
} else {
  cat("  (none)\n")
}
