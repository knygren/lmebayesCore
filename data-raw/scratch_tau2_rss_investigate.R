## Deep check: RSS_emp vs VarCorr-implied RSS — bug hunt + de-shrinkage
pkgload::load_all(".", quiet = TRUE)
pkgload::load_all("../lmebayes", quiet = TRUE)
library(lme4)
library(lme4) # nolint
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
k <- "distracted_ppvt"

lmer_fit <- lmer(form, data = dat, REML = TRUE)
vc <- extract_mer_variance_components(lmer_fit, design$re_coef_names)
tau2 <- vc$vcov_re[[k]]

lmer_fixef <- fixef(lmer_fit)
lmer_coef <- coef(lmer_fit)[["school_id"]]
lmer_ranef <- ranef(lmer_fit)[["school_id"]]

Xk <- design$X_hyper[[k]]
Zk <- design$Z[, k, drop = FALSE]
cat("X_hyper ppvt col names:", colnames(Xk), "\n")
cat("Z ppvt col name:", colnames(Zk), "\n\n")

## --- RSS_emp (current code path) ---
b_k <- lmer_coef[, k]
names(b_k) <- rownames(lmer_coef)
gk_lmer <- setNames(lmer_fixef["distracted_ppvt"], "(Intercept)")
b_al <- glmbayesCore:::.two_block_align_b_to_xhyper(b_k, Xk, gl)
eta <- as.numeric(Xk[gl, , drop = FALSE] %*% gk_lmer)
rss_emp <- sum((b_al - eta)^2)

## --- identities ---
cat("=== Identity checks (ppvt) ===\n")
cat("max |coef - (fixef + ranef)|:", max(abs(
  lmer_coef[, k] - lmer_fixef[k] - lmer_ranef[, k]
)), "\n")
cat("RSS_emp == sum(ranef^2):", all.equal(rss_emp, sum(lmer_ranef[, k]^2)), "\n")
cat("var(ranef):", var(lmer_ranef[, k]), "  VarCorr:", tau2, "\n")
cat("J * var(ranef):", J * var(lmer_ranef[, k]), "\n")
cat("J * VarCorr:", J * tau2, "\n\n")

## --- gamma from build_mu_all / case1 ICM ---
ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
fit1 <- lmerb(form, dat, pfamily_list(ps), ps$dispersion_ranef, simulate = FALSE)
g_icm <- fit1$fixef.mode[[k]]
b_icm <- fit1$ranef.mode[, k]
mu <- build_mu_all(design, fit1$fixef)$mu_all
mu_ppvt <- mu[k, ]
rss_icm <- sum((b_icm - mu_ppvt[gl])^2)
rss_icm2 <- sum((b_icm - as.numeric(Xk[gl, , drop = FALSE] %*% g_icm))^2)
cat("=== Using ICM / mu_all instead of lmer fixef ===\n")
cat("RSS (icm b, icm gamma):", rss_icm2, "\n")
cat("RSS (icm b, mu_all row):", rss_icm, "\n\n")

## --- lme4 internal random effects ---
u <- lme4::getME(lmer_fit, "u")  # unscaled?
b_lme4 <- lme4::getME(lmer_fit, "b")  # random effects in model scale
cat("=== lme4 getME ===\n")
cat("length(b):", length(b_lme4), "expected:", J * 3, "\n")
## map b to ppvt component — find index for distracted_ppvt in reTrms
re_idx <- lme4::findbars(form)[[1]]
cn <- colnames(lmer_ranef)
ppvt_col <- which(cn == k)
cat("ppvt column in ranef:", ppvt_col, "\n")
## b is stacked; use ranef as authoritative for ppvt per school
cat("max |getME b segment - ranef| (ppvt block): check via ranef only\n\n")

## --- observation-level: Z contribution to random slope variance ---
## For slope k: sum_i (z_ij * u_j)^2 contributions?
y <- dat$score_ppvt
groups <- dat$school_id
sigma2 <- vc$residual_var
## per-school n_j and sum of (distracted_ppvt * centered?) 
n_j <- as.numeric(table(groups)[gl])
names(n_j) <- gl

## REML-style: total SS of random slope effects weighted by design
## Var of BLUP = tau2 * (1 - lambda_j), lambda_j = sigma^2 / (sigma^2 + n_j * tau2 * v)
## For uncorrelated random slope with z = x_ij (distracted_ppvt value per obs):
x_ppvt <- dat$distracted_ppvt
## school-level: sum over obs in j of z_ij * u_j where u_j = ranef_j
u_j <- lmer_ranef[match(groups, gl), k]
contrib_obs <- u_j * x_ppvt  # not quite - Z matrix uses raw distracted_ppvt

## Shrinkage factor per school (random intercept/slope model approx for slope)
## lambda_j = 1 - var_cond / tau2 from lme4 condVar
re_cv <- ranef(lmer_fit, condVar = TRUE)[[1]]
if (!is.null(re_cv)) {
  cv_ppvt <- attr(re_cv, "postVar")
  if (!is.null(cv_ppvt)) {
    ## postVar is 3x3xJ for each school; extract ppvt diagonal
    idx <- match(k, cn)
    var_cond <- sapply(seq_len(J), function(j) {
      cv_ppvt[idx, idx, j]
    })
    names(var_cond) <- gl
    lambda <- 1 - var_cond / tau2
    cat("=== Shrinkage from condVar (ppvt) ===\n")
    cat("mean lambda:", mean(lambda), "  range:", range(lambda), "\n")
    cat("mean var_cond:", mean(var_cond), "  VarCorr:", tau2, "\n")
    ## de-shrunk ranef: u_unsh = ranef / sqrt(1-lambda)? or ranef / (1-lambda)?
    u_unsh <- lmer_ranef[, k] / sqrt(pmax(1 - lambda, 0.01))
    rss_unsh <- sum(u_unsh^2)
    cat("RSS if ranef / sqrt(1-lambda):", rss_unsh, "\n")
    cat("J * VarCorr:", J * tau2, "\n\n")
  }
}

## --- Alternative RSS: use coef totals around school-weighted mean ---
cat("=== Alternative RSS definitions ===\n")
cat("sum(coef^2):", sum(lmer_coef[, k]^2), "\n")
cat("sum((coef - mean(coef))^2):", sum((lmer_coef[, k] - mean(lmer_coef[, k]))^2), "\n")
cat("(J-1)*var(ranef):", (J - 1) * var(lmer_ranef[, k]), "\n")
cat("(J-1)*VarCorr:", (J - 1) * tau2, "\n")

## --- PT / trace: does lmer VarCorr satisfy J*tau2 = sum of something from fit? ---
cat("\n=== lme4 VarCorr vs sum of squares ===\n")
print(VarCorr(lmer_fit))
cat("\n")

## Manual: sigma and tau from lmer, average n per school
nbar <- mean(n_j)
lambda_bar <- sigma2 / (sigma2 + nbar * tau2)  # rough for scalar slope z=1
cat(sprintf("Rough lambda (z=1, nbar=%.1f): %.3f\n", nbar, lambda_bar))
cat(sprintf("var(BLUP) ≈ tau2*(1-lambda): %.2f * %.3f = %.2f (emp var ranef %.2f)\n",
            tau2, 1 - lambda_bar, tau2 * (1 - lambda_bar), var(lmer_ranef[, k])))
cat(sprintf("RSS_emp / (J*tau2) = %.3f\n", rss_emp / (J * tau2)))
cat(sprintf("If RSS_unshrink = RSS_emp / (1-lambda_bar): %.1f\n",
            rss_emp / (1 - lambda_bar)))
