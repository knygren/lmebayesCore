pkgload::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(lme4))
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
gl <- levels(design$groups)
J <- length(gl)
fit <- lmer(form, data = dat, REML = TRUE)
sigma2 <- summary(fit)$sigma^2
vc <- design$vcov_re
re <- design$re_coef_names
fixef_cal <- fixef_from_mer_fit(fit, design, "full_model")
b_cal <- ranef_from_mer_fit(fit, design, fixef_cal = fixef_cal)
mu_all <- build_mu_all(design, fixef_cal, gl)$mu_all
re_mat <- as.data.frame(ranef(fit)[["school_id"]])

cat("sigma^2 (residual):", round(sigma2, 4), "\n")
cat("var(y):", round(var(dat$score_ppvt), 4), "\n\n")

hdr <- sprintf(
  "%-18s %10s %10s %10s %10s %10s %10s %12s\n",
  "RE", "VarCorr", "tau2_rss", "var(b)", "var(u)",
  "VC/sigma2", "VC/tau2_rss", "VC*sigma2/tau2"
)
cat(hdr)
cat(strrep("-", nchar(trimws(hdr))), "\n", sep = "")

for (k in re) {
  tau_vc <- unname(vc[k])
  rss <- rss_k_from_mer_fit(b_cal[, k], mu_all[k, ], design$X_hyper[[k]], gl)
  tau_rss <- rss / (J - 1)
  vb <- var(b_cal[, k])
  vu <- var(re_mat[[k]])
  adj <- if (tau_rss > 0) tau_vc * sigma2 / tau_rss else NA_real_
  cat(sprintf(
    "%-18s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %12.4f\n",
    k, tau_vc, tau_rss, vb, vu, tau_vc / sigma2, tau_vc / tau_rss, adj
  ))
}

re_cv <- ranef(fit, condVar = TRUE)[["school_id"]]
pv <- attr(re_cv, "postVar")
cn <- colnames(re_mat)
cat("\n=== Adjusted VarCorr: VarCorr / data variance ===\n")
cat("  sigma2 (residual) = data variance in mixed-model sense\n")
cat("  var(y)            = total outcome variance\n\n")
cat(sprintf(
  "%-18s %10s %10s %10s %10s %10s %10s\n",
  "RE", "VarCorr", "tau2_rss", "VC/sigma2", "VC/var(y)",
  "tau2/sigma2", "VC/tau2_rss"
))
cat(strrep("-", 78), "\n", sep = "")
for (k in re) {
  tau_vc <- unname(vc[k])
  tau_rss <- rss_k_from_mer_fit(
    b_cal[, k], mu_all[k, ], design$X_hyper[[k]], gl
  ) / (J - 1)
  cat(sprintf(
    "%-18s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
    k, tau_vc, tau_rss,
    tau_vc / sigma2, tau_vc / var(dat$score_ppvt),
    tau_rss / sigma2, tau_vc / tau_rss
  ))
}

cat("\n=== Does adjusted VarCorr bridge to tau2_rss? ===\n")
cat("  adj_VC = VarCorr/sigma2  (ICC-like on residual scale)\n")
cat("  tau2_rss * (VarCorr/tau2_rss) = VarCorr  (trivial)\n")
cat("  tau2_rss * (VarCorr/var(BLUP)) should ≈ VarCorr when var(BLUP)=tau2_rss\n\n")
for (k in re) {
  tau_vc <- unname(vc[k])
  tau_rss <- rss_k_from_mer_fit(
    b_cal[, k], mu_all[k, ], design$X_hyper[[k]], gl
  ) / (J - 1)
  vb <- var(b_cal[, k])
  vu <- var(re_mat[[k]])
  cat(sprintf(
    "%-18s  VarCorr/var(BLUP)=%.2f  tau2_rss*ratio=%.2f  VarCorr=%.2f\n",
    k, tau_vc / vu, tau_rss * (tau_vc / vu), tau_vc
  ))
  cat(sprintf(
    "%-18s  VarCorr/var(b)=%.2f  tau2_rss*(VC/var(b))=%.2f\n",
    k, tau_vc / vb, tau_rss * (tau_vc / vb)
  ))
}

cat("\n=== De-shrink via condVar: tau2_rss / mean(1 - lambda) ===\n")
for (k in re) {
  idx <- match(k, cn)
  tau_vc <- unname(vc[k])
  var_cond <- vapply(seq_len(J), function(j) pv[[idx, idx, j]], numeric(1))
  lam <- 1 - var_cond / tau_vc
  tau_rss <- rss_k_from_mer_fit(
    b_cal[, k], mu_all[k, ], design$X_hyper[[k]], gl
  ) / (J - 1)
  deshrink <- tau_rss / mean(1 - lam)
  cat(sprintf(
    "%-18s mean_lambda=%.3f  deshrink=%.2f  VarCorr=%.2f  ratio=%.2f\n",
    k, mean(lam), deshrink, tau_vc, deshrink / tau_vc
  ))
}
