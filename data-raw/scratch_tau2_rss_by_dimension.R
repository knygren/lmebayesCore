## RSS per RE dimension: empirical (lmer b, gamma) vs VarCorr-implied
pkgload::load_all(".", quiet = TRUE)
pkgload::load_all("../lmebayes", quiet = TRUE)
library(lme4)
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
form <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

design <- model_setup(form, data = dat)
ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
pf_ing <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
re <- design$re_coef_names
gl <- levels(design$groups)
J <- length(gl)
n_obs <- nrow(dat)

lmer_fit <- lmer(form, data = dat, REML = TRUE)
vc <- extract_mer_variance_components(lmer_fit, re)
lmer_fixef <- fixef(lmer_fit)
lmer_coef <- coef(lmer_fit)[["school_id"]]
lmer_ranef <- ranef(lmer_fit)[["school_id"]]

## VarCorr from user table (matches lmer output)
tau2_lmer <- c(
  "(Intercept)"     = 199.29,
  distracted_ppvt = 22.54,
  distracted_a1   = 42.07
)

lmer_gamma_k <- function(k) {
  par <- colnames(design$X_hyper[[k]])
  if (k == "(Intercept)") {
    lmer_fixef[par]
  } else if (k == "distracted_ppvt") {
    stats::setNames(unname(lmer_fixef["distracted_ppvt"]), par)
  } else if (k == "distracted_a1") {
    stats::setNames(
      unname(c(
        lmer_fixef["distracted_a1"],
        lmer_fixef["free_reduced_lunch:distracted_a1"]
      )),
      par
    )
  }
}

rss_empirical <- function(k) {
  b_k <- lmer_coef[, k]
  names(b_k) <- rownames(lmer_coef)
  Xk <- design$X_hyper[[k]]
  gk <- lmer_gamma_k(k)
  b_al <- glmbayesCore:::.two_block_align_b_to_xhyper(b_k, Xk, gl)
  eta <- as.numeric(Xk[gl, , drop = FALSE] %*% gk)
  sum((b_al - eta)^2)
}

cat("Model: ", deparse(form), "\n", sep = "")
cat(sprintf("n_obs = %d, J = %d groups\n\n", n_obs, J))

cat(sprintf(
  "%-18s %10s %12s %12s %12s %12s %10s\n",
  "RE component", "VarCorr", "RSS_emp", "J*tau2", "(J-1)*tau2",
  "sum(ranef^2)", "RSS/VarCorr"
))
cat(strrep("-", 88), "\n", sep = "")

for (k in re) {
  tau2 <- tau2_lmer[k]
  rss <- rss_empirical(k)
  rss_ranef_only <- sum(lmer_ranef[, k]^2)
  cat(sprintf(
    "%-18s %10.2f %12.1f %12.1f %12.1f %12.1f %10.2f\n",
    k, tau2, rss, J * tau2, (J - 1) * tau2, rss_ranef_only, rss / tau2
  ))
}

cat("\nNotes:\n")
cat("  RSS_emp     = sum_j (b_j - X_j gamma)^2  with b = coef(lmer), gamma = fixef(lmer)\n")
cat("  J*tau2      = what RSS would be if J uncorrelated N(0,tau2) deviations (mean eta=0)\n")
cat("  (J-1)*tau2  = REML-style df for group variance\n")
cat("  sum(ranef^2) = sum_j ranef_j^2  (equals RSS_emp when X_k is intercept-only)\n")
cat("  RSS/VarCorr = effective df in RSS (empirical); ~J if unshrunk, <<J if shrunk BLUPs\n\n")

## ING closed-form tau2 at each RSS
cat("Closed-form tau2 from two_block_tau2_mode_ing (shape/rate from pf_ing):\n")
cat(sprintf("%-18s %10s %10s %10s\n", "RE", "RSS_emp", "tau2_lmer", "tau2_closed"))
cat(strrep("-", 52), "\n", sep = "")
for (k in re) {
  pl <- pf_ing[[k]]$prior_list
  b_k <- lmer_coef[, k]
  names(b_k) <- rownames(lmer_coef)
  closed <- two_block_tau2_mode_ing(
    b_k = b_k,
    X_k = design$X_hyper[[k]],
    gamma_k = lmer_gamma_k(k),
    shape = pl$shape,
    rate = pl$rate,
    group_levels = gl,
    method = "closed_form"
  )
  cat(sprintf(
    "%-18s %10.1f %10.2f %10.2f\n",
    k, rss_empirical(k), tau2_lmer[k], closed
  ))
}

cat("\nIf RSS were J*VarCorr instead of RSS_emp, closed tau2 would be:\n")
for (k in re) {
  pl <- pf_ing[[k]]$prior_list
  shape <- pl$shape[1L]
  rate <- pl$rate[1L]
  rss_vc <- J * tau2_lmer[k]
  tau2_cf <- (rate + rss_vc / 2) / (shape + J / 2 - 1)
  cat(sprintf("  %-18s RSS=%8.1f -> tau2_closed = %.2f (VarCorr = %.2f)\n",
              k, rss_vc, tau2_cf, tau2_lmer[k]))
}
