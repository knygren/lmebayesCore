## Scratch: timing for the bayesrules-based vignette fits.
suppressMessages(library(lmebayes))

cat("bayesrules installed:", requireNamespace("bayesrules", quietly = TRUE), "\n")
cat("glmmTMB installed:", requireNamespace("glmmTMB", quietly = TRUE), "\n")

data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) & !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c("score_ppvt", "distracted_a1", "distracted_ppvt",
                           "private_school", "title1", "free_reduced_lunch",
                           "school_id")])
)
cat("rows:", nrow(dat), " schools:", nlevels(droplevels(dat$school_id)), "\n")

form <- score_ppvt ~ private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 + free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

t0 <- system.time(ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01))
cat("Prior_Setup_GLMM:", t0[["elapsed"]], "s\n")

t1 <- system.time({
  fit <- lmerb(form, data = dat, pfamily_list = pfamily_list(ps),
               dispersion_ranef = ps$group.dispersion,
               n = 1000L, progbar = FALSE)
})
cat("lmerb n=1000:", t1[["elapsed"]], "s\n")
cat("object.size(fit):", format(object.size(fit), units = "MB"), "\n")
