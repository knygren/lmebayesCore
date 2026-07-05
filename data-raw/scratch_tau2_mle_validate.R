## MLE validation: ING tau^2 mode vs lmer VarCorr (big_word_club)
## Calls two_block_tau2_mode_ing() with lmer (b, gamma) and ING shape/rate only.
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

lmer_fit <- lme4::lmer(form, data = dat, REML = TRUE)
vc <- extract_mer_variance_components(lmer_fit, re)
lmer_fixef <- fixef(lmer_fit)
lmer_coef <- coef(lmer_fit)[["school_id"]]
tau2_plug <- two_block_tau2_plug_in_list(pf_ing, re)

## Map lmer fixef to Block~2 gamma_k per RE component (X_hyper column names).
lmer_gamma_k <- function(k, fe = lmer_fixef) {
  Xk <- design$X_hyper[[k]]
  par <- colnames(Xk)
  if (k == "(Intercept)") {
    fe[par]
  } else if (k == "distracted_ppvt") {
    stats::setNames(unname(fe["distracted_ppvt"]), par)
  } else if (k == "distracted_a1") {
    stats::setNames(
      unname(c(
        fe["distracted_a1"],
        fe["free_reduced_lunch:distracted_a1"]
      )),
      par
    )
  } else {
    stop("unknown RE component: ", k)
  }
}

cat("tau^2 plug-ins from pfamily_list:\n")
print(round(tau2_plug, 4))
cat("\ntau^2 validation at lmer inputs (coef + fixef -> gamma)\n\n")
cat(sprintf(
  "%-18s %10s %10s %10s %10s\n",
  "RE component", "lmer", "closed", "rglmb", "prior"
))
cat(strrep("-", 62), "\n", sep = "")
for (k in re) {
  Xk <- design$X_hyper[[k]]
  pl <- pf_ing[[k]]$prior_list
  b_k <- lmer_coef[, k]
  names(b_k) <- rownames(lmer_coef)
  gk <- lmer_gamma_k(k)
  closed <- two_block_tau2_mode_ing(
    b_k          = b_k,
    X_k          = Xk,
    gamma_k      = gk,
    shape        = pl$shape,
    rate         = pl$rate,
    group_levels = gl,
    method       = "closed_form"
  )
  rgl <- two_block_tau2_mode_ing(
    b_k          = b_k,
    X_k          = Xk,
    gamma_k      = gk,
    shape        = pl$shape,
    rate         = pl$rate,
    pfamily      = pf_ing[[k]],
    group_levels = gl,
    method       = "rglmb"
  )
  cat(sprintf(
    "%-18s %10.4f %10.4f %10.4f %10.4f\n",
    k, vc$vcov_re[[k]], closed, rgl$tau2, tau2_plug[[k]]
  ))
}

cat("\nBlock~1 conditional prior Sigma_ranef at lmer tau^2 plug-ins:\n")
Sigma_lmer <- two_block_Sigma_ranef(vc$vcov_re, re)
print(round(Sigma_lmer, 4))

gamma_lmer <- lapply(re, lmer_gamma_k)
names(gamma_lmer) <- re
prior_b <- two_block_conditional_prior_ranef(
  design, gamma_lmer, tau2 = vc$vcov_re
)
cat("\nConditional prior mu_all (first 3 schools, intercept RE):\n")
print(round(prior_b$mu_all[1, 1:3], 4))

post_b <- two_block_conditional_posterior_ranef(
  design, gamma_lmer, tau2 = vc$vcov_re,
  sigma2 = vc$residual_var, group_level = gl[1]
)
cat("\nPosterior b for first school (mean, diag Sigma):\n")
print(round(post_b$mean, 4))
print(round(diag(post_b$Sigma), 4))
