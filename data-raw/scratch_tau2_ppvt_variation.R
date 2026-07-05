## Diagnose tau2: lmer VarCorr vs BLUP variation vs closed-form RSS
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
gl <- levels(design$groups)
J <- length(gl)

lmer_fit <- lmer(form, data = dat, REML = TRUE)
vc_df <- as.data.frame(VarCorr(lmer_fit))
vc_re <- vc_df$vcov[vc_df$grp == "school_id"]
names(vc_re) <- vc_df$var1[vc_df$grp == "school_id"]

lmer_fixef <- fixef(lmer_fit)
lmer_coef <- coef(lmer_fit)[["school_id"]]
lmer_ranef <- ranef(lmer_fit)[["school_id"]]

lmer_gamma_k <- function(k) {
  Xk <- design$X_hyper[[k]]
  par <- colnames(Xk)
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

tau2_from_rss <- function(rss, shape, rate) {
  (rate + rss / 2) / (shape + J / 2 - 1)
}

cat("Compare lmer VarCorr to BLUP sample variance and closed-form tau2\n")
cat(sprintf(
  "%-18s %8s %8s %8s %8s %8s\n",
  "RE", "VarCorr", "var(ranef)", "sd(coef)", "RSS", "closed"
))
cat(strrep("-", 66), "\n", sep = "")
for (k in design$re_coef_names) {
  b_tot <- lmer_coef[, k]
  names(b_tot) <- rownames(lmer_coef)
  b_ran <- lmer_ranef[, k]
  gk <- lmer_gamma_k(k)
  Xk <- design$X_hyper[[k]]
  pl <- pf_ing[[k]]$prior_list
  rss <- sum((
    b_tot - as.numeric(Xk[gl, , drop = FALSE] %*% gk)
  )^2)
  closed <- two_block_tau2_mode_ing(
    b_k = b_tot, X_k = Xk, gamma_k = gk,
    shape = pl$shape, rate = pl$rate,
    group_levels = gl, method = "closed_form"
  )
  cat(sprintf(
    "%-18s %8.2f %8.2f %8.2f %8.1f %8.2f\n",
    k, vc_re[k], var(b_ran), stats::sd(b_tot), rss, closed
  ))
}

k <- "distracted_ppvt"
cat("\n--- ppvt: totals look wide, ranefs do not ---\n")
b_tot <- lmer_coef[, k]
b_ran <- lmer_ranef[, k]
cat(sprintf(
  "fixef ppvt = %.3f; coef range [%.2f, %.2f]; ranef range [%.2f, %.2f]\n",
  lmer_fixef["distracted_ppvt"], min(b_tot), max(b_tot), min(b_ran), max(b_ran)
))
cat(sprintf(
  "VarCorr / var(ranef) = %.2f / %.2f = %.1fx\n",
  vc_re[k], var(b_ran), vc_re[k] / var(b_ran)
))
cat(
  "\nClosed-form tau2 uses RSS = sum((b - X gamma)^2) = sum(ranef^2).\n",
  "That matches var(BLUP) ~ 4.3, not lmer VarCorr ~ 22.5 (latent RE variance).\n",
  "Treating shrunk BLUPs as direct N(eta, tau2) observations is the mismatch.\n",
  sep = ""
)

