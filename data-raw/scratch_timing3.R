## Scratch: timing for the GLMM (Poisson / binomial) vignette fits.
suppressMessages(library(lmebayes))

## --- Poisson: airbnb ------------------------------------------------------
data(airbnb, package = "bayesrules")
ab <- airbnb[, c("neighborhood", "reviews", "rating", "room_type")]
ab <- ab[complete.cases(ab), ]
ab$neighborhood <- factor(ab$neighborhood)
ab$rating_c <- ab$rating - mean(ab$rating)
cat("airbnb rows:", nrow(ab), " neighborhoods:", nlevels(ab$neighborhood), "\n")

form_p <- reviews ~ rating_c + (1 + rating_c || neighborhood)
t0 <- system.time(
  ps_p <- Prior_Setup_GLMM(form_p, data = ab, family = poisson(),
                           pop.pwt = 0.01)
)
cat("Prior_Setup_GLMM poisson:", t0[["elapsed"]], "s\n")

t1 <- system.time({
  fit_p <- glmerb(form_p, data = ab, family = poisson(),
                  pfamily_list = pfamily_list(ps_p), n = 100L,
                  progbar = FALSE)
})
cat("glmerb poisson n=100:", t1[["elapsed"]], "s\n")
cat("m_convergence:", fit_p$m_convergence, "\n")
cat("size:", format(object.size(fit_p), units = "MB"), "\n")

## --- Binomial: book_banning ----------------------------------------------
data(book_banning, package = "bayesrules")
bb <- book_banning[, c("state", "removed", "violent")]
bb <- bb[complete.cases(bb), ]
bb$state <- factor(bb$state)
cat("book_banning rows:", nrow(bb), " states:", nlevels(bb$state), "\n")

form_b <- removed ~ violent + (1 + violent || state)
t2 <- system.time(
  ps_b <- Prior_Setup_GLMM(form_b, data = bb, family = binomial(),
                           pop.pwt = 0.01)
)
cat("Prior_Setup_GLMM binomial:", t2[["elapsed"]], "s\n")

t3 <- system.time({
  fit_b <- glmerb(form_b, data = bb, family = binomial(),
                  pfamily_list = pfamily_list(ps_b), n = 100L,
                  progbar = FALSE)
})
cat("glmerb binomial n=100:", t3[["elapsed"]], "s\n")
cat("m_convergence:", fit_b$m_convergence, "\n")
cat("size:", format(object.size(fit_b), units = "MB"), "\n")
