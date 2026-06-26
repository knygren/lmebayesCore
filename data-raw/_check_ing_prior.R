library(lmebayes)
library(bayesrules)
data(big_word_club)
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) &
    invalid_ppvt == 0L
)
form <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

ps <- Prior_Setup_lmebayes(
  form,
  data           = dat,
  pwt            = 0.01,
  pwt_dispersion = 0.2
)
pf <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
J <- nlevels(ps$design$groups)
cat("J =", J, "\n\n")
for (k in names(pf)) {
  pl <- pf[[k]]$prior_list
  p_k <- length(pl$mu)
  cat(k, " (p_k =", p_k, ")\n")
  cat("  shape =", pl$shape, " rate =", pl$rate, "\n")
  cat("  dispersion_fixef =", ps$prior_list[[k]]$dispersion_fixef, "\n")
  cat("  prior mean tau2 = rate/(shape-1) =", pl$rate / (pl$shape - 1), "\n")
  cat("  disp_lower =", pl$disp_lower, " disp_upper =", pl$disp_upper, "\n")
  cat("  n_prior_disp =", ps$n_prior_dispersion[[k]], "\n")
  cat("  implied n_prior = 2*shape - 1 - p =", 2 * pl$shape - 1 - p_k, "\n\n")
}

## One Block-2 hyper-regression draw (school random intercepts as y)
design <- ps$design
k <- "(Intercept)"
X_k <- design$X_hyper[[k]]
set.seed(1)
b_fake <- rnorm(nrow(X_k), mean = 100, sd = 10)
cat("=== test rglmb on", k, "===\n")
tryCatch(
  {
    fit <- glmbayesCore::rglmb(
      n = 1L,
      y = b_fake,
      x = X_k,
      family = gaussian(),
      pfamily = pf[[k]],
      verbose = TRUE,
      use_parallel = FALSE
    )
    cat("OK dispersion =", fit$dispersion, "\n")
  },
  error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  }
)
