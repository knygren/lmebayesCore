## Prior_Setup(shape_df) demo: matching lm()-style variance (n - p) vs n-counting
##
## Prior means = full-model MLE (intercept_source / effects_source = "full_model")
## so comparisons focus on variance, not location.  n = 20 PlantGrowth-style rows.
##
## Idea (weak prior on coefficients, see ?Prior_Setup and shape_df):
##   * lm()-like (residual df n - p): use dNormal_Gamma with default shape_df
##     (shape ~ n_prior/2) together with dIndependent_Normal_Gamma and shape_df
##     "n_prior+p" (and the same rate rule as in the package).
##   * Variance analogous to counting n (not n - p): use dIndependent_Normal_Gamma
##     with default shape_df together with dNormal_Gamma and shape_df "n_prior-p".
##
## Calibrations used here:
##   * For the n-p style (section 1): pwt = 0.01 for both Prior_Setup calls; do not
##     pass n_prior (it is implied by pwt).  First call: default shape_df for NG;
##     second: shape_df = "n_prior+p" for ING.
##     pwt = 0.01 is not a flat prior: implied n_prior = (pwt/(1-pwt))*n_effective
##     (about 0.2 when n = 20), plus a Gamma prior on residual precision.  Posterior
##     variance can sit slightly *below* vcov(lm) (ratio < 1 on the diagonal), not
##     exactly equal; that is expected unless pwt -> 0.
##   * For the n style (section 2): n_prior_sec2 = 3 for both Prior_Setup calls
##     (requires n_prior > p for "n_prior-p").  Section 2 replicates each row 5 times
##     (n_effective = 100).  Gamma *rates* for lmb: ING has shape = n_prior/2 = 1.5,
##     so rate = (shape - 1) * d_np with d_np = summary(lm)$sigma^2 (RSS/(n-p)).
##     NG with shape_df "n_prior-p" has shape = (3-2)/2 = 0.5 < 1, so E[sigma^2]
##     calibration does not apply; section 2 uses Prior_Setup default rate for NG.
##   * Sections 3-4 mirror 1-2 with disp_type = "Post_mean": Nelder-Mead adjusts
##     (dispersion, lambda) so the *conjugate fragment's* identity for E[dispersion|y]
##     holds jointly with the (1-pwt)*MLE + pwt*mu coefficient target.  That pins
##     Prior_Setup$dispersion (and Sigma) used in lmb(); it does *not* force
##     vcov(lmb) to equal vcov(lm) at finite pwt.
##     Typical pattern (full_model mu, this script): Sec 3 NG often leaves d* equal to
##     starting d0 = RSS/(n-2) which matches summary(lm)$sigma^2 here (p=2); Sec 3 ING
##     lowers d* vs lm sigma^2 (larger Gamma shape from n_prior+p).  Sec 4 ING often
##     gives d* = summary(lm)$sigma^2 on replicated data; Sec 4 NG can nudge d* slightly
##     above d_np.  Ratios vcov(lmb)/vcov(lm) often stay below 1 on the diagonal; the
##     (1-pwt)*vcov(lm) heuristic still tracks Sec 1 NG and Sec 2 ING well; Post_mean
##     rows are not uniformly closer to raw lm than OLS_mean.
##
## Why Post_mean can move mean-diagonal vcov(lmb)/vcov(lm) *farther* from 1 than OLS_mean:
##   * Nelder-Mead minimizes mismatch to the *conjugate fragment* (fixed-point dispersion
##     and (1-pwt)*MLE+pwt*mu for beta in that fragment), not mismatch of *marginal*
##     posterior Var(beta) from lmb() to vcov(lm).
##   * Post_mean rewrites Sigma (lambda) and dispersion together; lmb() then integrates
##     over tau (and uses the full NG / ING structure). A good d* for the fragment need
##     not imply the same marginal beta spread as the OLS_mean prior parameterization.
##   * Example: Sec 3 ING lowers d* vs summary(lm)$sigma^2, which pulls uncertainty down
##     vs the OLS_mean path with the same pwt. Sec 3 NG often leaves d* at the OLS_mean
##     plug-in but still changes nothing if lambda is flat; small MC noise can reorder
##     ratios slightly vs Sec 1.
##
## End of script: summary table - vcov ratios, target_df_denom, post_mean_disp,
## disp_vs_lm = E_post[sigma^2]/dispersion_lm.  post_mean_prec = E_post[1/sigma^2] via
## mean(1/dispersion draws).  prec_vs_lm = post_mean_prec / (1/dispersion_lm); not 1/disp_vs_lm.
## Summary table includes n_over_n_plus_nprior = n/(n+n_prior) between vs_lm and vs_scaled,
## then prior_setup_dispersion = Prior_Setup()$dispersion (Gaussian plug-in / Post_mean optimum).
##
## MC draws: n_mc = 100000 for stable vcov(lmb).  Early output includes a "scaled
## vcov(lm)" = vcov(lm) * (1 - pwt): heuristic when posterior beta precision scales
## like 1/(1-pwt) under the Zellner setup (Sections 1,3: pwt = 0.01; Sections 2,4: pwt
## from n_prior_sec2 and n_effective = 100 after 5-fold replication).  Congruence summaries use
## M = U^{-T} vcov(lmb) U^{-1} with vcov(lm) = U'U (chol); eigenvalues of M replace
## elementwise ratios for a matrix-consistent scale comparison.
##
## Run: demo(Ex_10_Prior_Setup_shape_df, package = "glmbayes")

library(glmbayes)

## Matrix comparison of covariances: V_base = crossprod(U) with U = chol(V_base)
## upper-triangular.  M = U^{-T} V_post U^{-1} = t(Ui) %*% V_post %*% Ui with
## Ui = U^{-1}.  Eigenvalues of M are generalized eigenvalues of (V_post, V_base);
## if V_post = k * V_base then M = k * I (all eigenvalues k).
vcov_congruence <- function(V_post, V_base) {
  U <- chol(V_base)
  Ui <- backsolve(U, diag(nrow(U)))
  M <- crossprod(Ui, V_post %*% Ui)
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  list(M = M, eigenvalues = ev)
}

congruence_line <- function(V_post, label, V_base, base_label) {
  cq <- vcov_congruence(V_post, V_base)
  ev <- cq$eigenvalues
  sprintf(
    "%s vs %s: eigenvalues min / mean / max = %s / %s / %s",
    label,
    base_label,
    format(min(ev), digits = 6),
    format(mean(ev), digits = 6),
    format(max(ev), digits = 6)
  )
}

## Posterior mean dispersion (sigma^2) from lmb: mean of MC draws in $dispersion.
post_mean_lmb_dispersion <- function(fit) {
  d <- fit$dispersion
  if (length(d) == 1L) as.numeric(d) else mean(as.numeric(d))
}

## Posterior expected data precision E[tau | y], tau = 1/sigma^2: MC mean of 1/dispersion.
## Not 1/E[sigma^2] and not derived from disp_vs_lm (Jensen: E[1/d] != 1/E[d]).
post_mean_lmb_precision <- function(fit) {
  d <- fit$dispersion
  if (length(d) == 1L) {
    d1 <- as.numeric(d)
    if (!is.finite(d1) || d1 <= 0) stop("lmb dispersion must be finite and positive.")
    return(1 / d1)
  }
  dd <- as.numeric(d)
  if (any(!is.finite(dd) | dd <= 0)) stop("lmb dispersion draws must be finite and positive.")
  mean(1 / dd)
}

ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)

## Section 2 only: replicate each observation 5 times (n = 100) for smaller implied pwt.
rep_sec2 <- 5L
idx_s2 <- rep(seq_along(weight), each = rep_sec2)
weight_s2 <- weight[idx_s2]
group_s2 <- group[idx_s2]
n_eff_s2 <- length(weight_s2)
n_prior_sec2 <- 3L

n_mc <- 100000L

fit_lm <- lm(weight ~ group, x = TRUE, y = TRUE)
V_lm <- vcov(fit_lm)

fit_lm_s2 <- lm(weight_s2 ~ group_s2, x = TRUE, y = TRUE)
V_lm_s2 <- vcov(fit_lm_s2)

pwt_sec1 <- 0.01
V_lm_scaled_s1 <- V_lm * (1 - pwt_sec1)
pwt_sec2 <- n_prior_sec2 / (n_prior_sec2 + n_eff_s2)
V_lm_scaled_s2 <- V_lm_s2 * (1 - pwt_sec2)

cat("\n======== Reference: vcov(lm) and scaled vcov (Zellner variance heuristic) =========\n")
cat(
  "Heuristic: if beta precision scales by 1/(1-pwt), then Var ~ (1-pwt)*vcov(lm).\n",
  "So scaled vcov(lm) = vcov(lm) * (1-pwt).  Compare posterior vcov(lmb) to these.\n\n",
  sep = ""
)
cat("--- Ordinary least squares ---\n")
print(stats::coef(fit_lm))
cat("\nvcov(lm):\n")
print(V_lm)
cat(
  "\n--- Section 1 (pwt = ", pwt_sec1, "): scaled = vcov(lm) * (1-pwt) = vcov(lm) * ",
  format(1 - pwt_sec1, digits = 6), " ---\n",
  sep = ""
)
print(V_lm_scaled_s1)
cat(
  "\n--- Section 2 (5x replication: n_effective = ", n_eff_s2,
  "; n_prior = ", n_prior_sec2, " => pwt = ", format(pwt_sec2, digits = 6),
  "): lm on replicated data; scaled = vcov(lm_s2) * (1-pwt) = vcov(lm_s2) * ",
  format(1 - pwt_sec2, digits = 6), " ---\n",
  sep = ""
)
cat("\nvcov(lm) on replicated data (same coef as n=20, smaller SE^2):\n")
print(V_lm_s2)
print(V_lm_scaled_s2)

## ----- Section 1: lm-like (n - p) variance -----
## Prior_Setup twice: pwt = 0.01, no n_prior; NG uses default shape_df, ING uses n_prior+p.

cat("\n")
cat("######################################################################\n")
cat("# Section 1: priors aimed at lm()-style (n - p) residual counting   #\n")
cat("#   dNormal_Gamma: default shape_df (n_prior)                       #\n")
cat("#   dIndependent_Normal_Gamma: shape_df = n_prior+p, pwt = 0.01    #\n")
cat("######################################################################\n")

p1_ng <- Prior_Setup(
  weight ~ group,
  disp_type = "OLS_mean",
  pwt = 0.01,
  shape_df = "n_prior",
  intercept_source = "full_model",
  effects_source = "full_model"
)

p1_ing <- Prior_Setup(
  weight ~ group,
  disp_type = "OLS_mean",
  pwt = 0.01,
  shape_df = "n_prior+p",
  intercept_source = "full_model",
  effects_source = "full_model"
)

fit_p1_ng <- lmb(
  weight ~ group,
  pfamily = dNormal_Gamma(
    p1_ng$mu,
    p1_ng$Sigma / p1_ng$dispersion,
    shape = p1_ng$shape,
    rate = p1_ng$rate
  ),
  n = n_mc
)

fit_p1_ing <- lmb(
  weight ~ group,
  pfamily = dIndependent_Normal_Gamma(
    p1_ing$mu,
    p1_ing$Sigma,
    shape = p1_ing$shape,
    rate = p1_ing$rate
  ),
  n = n_mc
)

cat("\n======== Section 1: location (lm vs prior means) =========\n")
print(stats::coef(fit_lm))
print(cbind(lm = coef(fit_lm), ng_mu = p1_ng$mu[, 1], ing_mu = p1_ing$mu[, 1]))

cat("\n======== Section 1: posterior vcov(lmb), dNormal_Gamma (default shape_df) =========\n")
print(vcov(fit_p1_ng))

cat("\n======== Section 1: posterior vcov(lmb), dIndependent_Normal_Gamma (n_prior+p) =========\n")
print(vcov(fit_p1_ing))

cat("\n======== Section 1: ratios vcov(lmb) / vcov(lm) =========\n")
R1_ng <- vcov(fit_p1_ng) / V_lm
R1_ing <- vcov(fit_p1_ing) / V_lm
cat("dNormal_Gamma (default n_prior):\n")
print(R1_ng)
cat("dIndependent_Normal_Gamma (n_prior+p):\n")
print(R1_ing)
cat(sprintf(
  "Mean diagonal ratio vs lm:  NG (default) %5.3f;  ING (n_prior+p) %5.3f\n",
  mean(diag(R1_ng)),
  mean(diag(R1_ing))
))
R1_ng_s <- vcov(fit_p1_ng) / V_lm_scaled_s1
R1_ing_s <- vcov(fit_p1_ing) / V_lm_scaled_s1
cat("\nSection 1: vcov(lmb) / scaled vcov(lm)  (near 1 if Zellner scalar heuristic fits):\n")
cat("dNormal_Gamma (default n_prior):\n")
print(R1_ng_s)
cat("dIndependent_Normal_Gamma (n_prior+p):\n")
print(R1_ing_s)
cat(sprintf(
  "Mean diagonal vs scaled:  NG %5.3f;  ING (n_prior+p) %5.3f\n",
  mean(diag(R1_ng_s)),
  mean(diag(R1_ing_s))
))
cat(
  "\nNote: ratios vs raw lm slightly below 1 are normal (finite pwt, prior on tau). ",
  "ING + n_prior+p often deviates more from the pure (1-pwt) variance story than NG.\n",
  sep = ""
)

cat("\nSection 1: congruence M = U^{-T} vcov(lmb) U^{-1}  (V_lm = U'U from chol):\n")
cat(congruence_line(vcov(fit_p1_ng), "NG", V_lm, "vcov(lm)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p1_ing), "ING (n_prior+p)", V_lm, "vcov(lm)"), "\n", sep = "")
cat(
  "  (If V_post = k * V_lm, all eigenvalues equal k; vs scaled target k = (1-pwt), compare to NG.)\n",
  sep = ""
)
cat(congruence_line(vcov(fit_p1_ng), "NG", V_lm_scaled_s1, "scaled vcov(lm)"), "\n", sep = "")
cat(
  congruence_line(vcov(fit_p1_ing), "ING (n_prior+p)", V_lm_scaled_s1, "scaled vcov(lm)"),
  "\n",
  sep = ""
)

## ----- Section 2: n-counting style variance -----
## Prior_Setup twice: n_prior_sec2; ING default shape_df, NG n_prior-p.
## Gamma rates for lmb: ING E[sigma^2] prior = rate/(shape-1) with d_np = RSS/(n-p);
## NG uses default rate when shape <= 1.

cat("\n")
cat("######################################################################\n")
cat("# Section 2: priors aimed at n-style (not n - p) counting          #\n")
cat("#   Data: each plant weight repeated ", rep_sec2, " times (n = ", n_eff_s2, ")   #\n", sep = "")
cat("#   n_prior = ", n_prior_sec2, "; ING: default shape_df; NG: n_prior-p           #\n", sep = "")
cat("#   ING: Gamma rate = (shape - 1) * d_np; NG: shape < 1 => default Prior_Setup rate. #\n")
cat("######################################################################\n")

p2_ing <- Prior_Setup(
  weight_s2 ~ group_s2,
  disp_type = "OLS_mean",
  n_prior = n_prior_sec2,
  shape_df = "n_prior",
  intercept_source = "full_model",
  effects_source = "full_model"
)

p2_ng <- Prior_Setup(
  weight_s2 ~ group_s2,
  disp_type = "OLS_mean",
  n_prior = n_prior_sec2,
  shape_df = "n_prior-p",
  intercept_source = "full_model",
  effects_source = "full_model"
)

## Prior on precision tau ~ Gamma(shape, rate) (glmbayes parameterization).
## sigma^2 = 1/tau has prior mean E[sigma^2] = rate / (shape - 1) for shape > 1.
## Classical unbiased residual variance: d_np = RSS/(n-p) = summary(lm)$sigma^2.
## Choose rate so E[sigma^2]_prior = d_np  =>  rate = (shape - 1) * d_np.
## Different shape_df => different shape => different rate (Prior_Setup $rate ignored here).
d_np_s2 <- as.numeric(summary(fit_lm_s2)$sigma^2)
shape_p2_ing <- p2_ing$shape
shape_p2_ng <- p2_ng$shape
rate_p2_ing <- (shape_p2_ing - 1) * d_np_s2
rate_p2_ng <- if (shape_p2_ng > 1) {
  (shape_p2_ng - 1) * d_np_s2
} else {
  p2_ng$rate
}

cat("\n======== Section 2: Gamma rate calibration (ING); NG default if shape <= 1 =========\n")
cat("d_np = summary(lm_s2)$sigma^2 (RSS/(n-p)) =", format(d_np_s2, digits = 8), "\n")
cat("ING: shape =", format(shape_p2_ing, digits = 6), "  rate = (shape-1)*d_np =", format(rate_p2_ing, digits = 8), "\n")
cat("  (Prior_Setup rate was", format(p2_ing$rate, digits = 8), ")\n")
if (shape_p2_ng > 1) {
  cat("NG:  shape =", format(shape_p2_ng, digits = 6), "  rate = (shape-1)*d_np =", format(rate_p2_ng, digits = 8), "\n")
} else {
  cat(
    "NG:  shape =", format(shape_p2_ng, digits = 6),
    "  (<= 1: no E[sigma^2] calibration)  rate = Prior_Setup default =", format(rate_p2_ng, digits = 8), "\n"
  )
}
cat("  (Prior_Setup rate was", format(p2_ng$rate, digits = 8), ")\n")

fit_p2_ing <- lmb(
  weight_s2 ~ group_s2,
  pfamily = dIndependent_Normal_Gamma(
    p2_ing$mu,
    p2_ing$Sigma,
    shape = shape_p2_ing,
    rate = rate_p2_ing
  ),
  n = n_mc
)

fit_p2_ng <- lmb(
  weight_s2 ~ group_s2,
  pfamily = dNormal_Gamma(
    p2_ng$mu,
    p2_ng$Sigma / p2_ng$dispersion,
    shape = shape_p2_ng,
    rate = rate_p2_ng
  ),
  n = n_mc
)

cat("\n======== Section 2: location (lm vs prior means, replicated n) =========\n")
print(stats::coef(fit_lm_s2))
print(cbind(lm = coef(fit_lm_s2), ing_mu = p2_ing$mu[, 1], ng_mu = p2_ng$mu[, 1]))

cat("\n======== Section 2: reference vcov(lm) on replicated data =========\n")
print(V_lm_s2)

cat("\n======== Section 2: posterior vcov(lmb), dIndependent_Normal_Gamma (default shape_df) =========\n")
print(vcov(fit_p2_ing))

cat("\n======== Section 2: posterior vcov(lmb), dNormal_Gamma (n_prior-p) =========\n")
print(vcov(fit_p2_ng))

cat("\n======== Section 2: ratios vcov(lmb) / vcov(lm) on replicated data =========\n")
R2_ing <- vcov(fit_p2_ing) / V_lm_s2
R2_ng <- vcov(fit_p2_ng) / V_lm_s2
cat("dIndependent_Normal_Gamma (default shape_df, calibrated rate):\n")
print(R2_ing)
cat("dNormal_Gamma (n_prior-p; Prior_Setup default rate, shape = 0.5):\n")
print(R2_ng)
cat(sprintf(
  "Mean diagonal ratio vs lm:  ING (default) %5.3f;  NG (n_prior-p) %5.3f\n",
  mean(diag(R2_ing)),
  mean(diag(R2_ng))
))
R2_ing_s <- vcov(fit_p2_ing) / V_lm_scaled_s2
R2_ng_s <- vcov(fit_p2_ng) / V_lm_scaled_s2
cat("\nSection 2: vcov(lmb) / scaled vcov(lm_s2)  (pwt = n_prior/(n_prior+n_effective)):\n")
cat("dIndependent_Normal_Gamma (default n_prior):\n")
print(R2_ing_s)
cat("dNormal_Gamma (n_prior-p; default rate):\n")
print(R2_ng_s)
cat(sprintf(
  "Mean diagonal vs scaled:  ING %5.3f;  NG (n_prior-p) %5.3f\n",
  mean(diag(R2_ing_s)),
  mean(diag(R2_ng_s))
))

cat("\nSection 2: congruence M = U^{-T} vcov(lmb) U^{-1}  (base = vcov(lm_s2)):\n")
cat(congruence_line(vcov(fit_p2_ing), "ING (default)", V_lm_s2, "vcov(lm_s2)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p2_ng), "NG (n_prior-p)", V_lm_s2, "vcov(lm_s2)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p2_ing), "ING (default)", V_lm_scaled_s2, "scaled vcov(lm_s2)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p2_ng), "NG (n_prior-p)", V_lm_scaled_s2, "scaled vcov(lm_s2)"), "\n", sep = "")

## ----- Section 3: same as section 1 but disp_type = "Post_mean" -----
## Prior_Setup runs Nelder-Mead so E[dispersion|y] matches the conjugate fragment;
## returned dispersion, Sigma, rate are used as-is in lmb().

cat("\n")
cat("######################################################################\n")
cat("# Section 3: like Section 1, disp_type = Post_mean (fixed-point d)   #\n")
cat("#   dNormal_Gamma: default shape_df; ING: shape_df = n_prior+p       #\n")
cat("#   pwt = 0.01; compare posterior dispersion / vcov to classical lm  #\n")
cat("######################################################################\n")

p3_ng <- Prior_Setup(
  weight ~ group,
  disp_type = "Post_mean",
  pwt = 0.01,
  shape_df = "n_prior",
  intercept_source = "full_model",
  effects_source = "full_model"
)

p3_ing <- Prior_Setup(
  weight ~ group,
  disp_type = "Post_mean",
  pwt = 0.01,
  shape_df = "n_prior+p",
  intercept_source = "full_model",
  effects_source = "full_model"
)

fit_p3_ng <- lmb(
  weight ~ group,
  pfamily = dNormal_Gamma(
    p3_ng$mu,
    p3_ng$Sigma / p3_ng$dispersion,
    shape = p3_ng$shape,
    rate = p3_ng$rate
  ),
  n = n_mc
)

fit_p3_ing <- lmb(
  weight ~ group,
  pfamily = dIndependent_Normal_Gamma(
    p3_ing$mu,
    p3_ing$Sigma,
    shape = p3_ing$shape,
    rate = p3_ing$rate
  ),
  n = n_mc
)

cat("\n======== Section 3: Prior_Setup dispersion (Post_mean fixed point) vs lm =========\n")
cat(
  "summary(lm)$sigma^2 (RSS/(n-p)) =", format(as.numeric(summary(fit_lm)$sigma^2), digits = 8), "\n"
)
cat("dNormal_Gamma:        Prior_Setup$dispersion =", format(p3_ng$dispersion, digits = 8), "\n")
cat("dIndep_Normal_Gamma:  Prior_Setup$dispersion =", format(p3_ing$dispersion, digits = 8), "\n")

cat("\n======== Section 3: posterior vcov(lmb) =========\n")
cat("dNormal_Gamma (default shape_df, Post_mean):\n")
print(vcov(fit_p3_ng))
cat("\ndIndependent_Normal_Gamma (n_prior+p, Post_mean):\n")
print(vcov(fit_p3_ing))

cat("\n======== Section 3: ratios vcov(lmb) / vcov(lm) =========\n")
R3_ng <- vcov(fit_p3_ng) / V_lm
R3_ing <- vcov(fit_p3_ing) / V_lm
cat("dNormal_Gamma:\n")
print(R3_ng)
cat("dIndependent_Normal_Gamma (n_prior+p):\n")
print(R3_ing)
cat(sprintf(
  "Mean diagonal ratio vs lm:  NG %5.3f;  ING (n_prior+p) %5.3f\n",
  mean(diag(R3_ng)),
  mean(diag(R3_ing))
))
R3_ng_s <- vcov(fit_p3_ng) / V_lm_scaled_s1
R3_ing_s <- vcov(fit_p3_ing) / V_lm_scaled_s1
cat("\nSection 3: vcov(lmb) / scaled vcov(lm):\n")
cat(sprintf(
  "Mean diagonal vs scaled:  NG %5.3f;  ING (n_prior+p) %5.3f\n",
  mean(diag(R3_ng_s)),
  mean(diag(R3_ing_s))
))

cat("\nSection 3: congruence vs vcov(lm) / scaled:\n")
cat(congruence_line(vcov(fit_p3_ng), "NG Post_mean", V_lm, "vcov(lm)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p3_ing), "ING Post_mean", V_lm, "vcov(lm)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p3_ng), "NG Post_mean", V_lm_scaled_s1, "scaled vcov(lm)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p3_ing), "ING Post_mean", V_lm_scaled_s1, "scaled vcov(lm)"), "\n", sep = "")

## ----- Section 4: same as section 2 but disp_type = "Post_mean" -----
## Use Prior_Setup outputs (no manual Gamma rate override).

cat("\n")
cat("######################################################################\n")
cat("# Section 4: like Section 2, disp_type = Post_mean                   #\n")
cat("#   Replicated n = 100, n_prior = 3; ING default / NG n_prior-p      #\n")
cat("######################################################################\n")

p4_ing <- Prior_Setup(
  weight_s2 ~ group_s2,
  disp_type = "Post_mean",
  n_prior = n_prior_sec2,
  shape_df = "n_prior",
  intercept_source = "full_model",
  effects_source = "full_model"
)

p4_ng <- Prior_Setup(
  weight_s2 ~ group_s2,
  disp_type = "Post_mean",
  n_prior = n_prior_sec2,
  shape_df = "n_prior-p",
  intercept_source = "full_model",
  effects_source = "full_model"
)

fit_p4_ing <- lmb(
  weight_s2 ~ group_s2,
  pfamily = dIndependent_Normal_Gamma(
    p4_ing$mu,
    p4_ing$Sigma,
    shape = p4_ing$shape,
    rate = p4_ing$rate
  ),
  n = n_mc
)

fit_p4_ng <- lmb(
  weight_s2 ~ group_s2,
  pfamily = dNormal_Gamma(
    p4_ng$mu,
    p4_ng$Sigma / p4_ng$dispersion,
    shape = p4_ng$shape,
    rate = p4_ng$rate
  ),
  n = n_mc
)

cat("\n======== Section 4: Prior_Setup dispersion (Post_mean) vs lm_s2 =========\n")
cat(
  "summary(lm_s2)$sigma^2 =", format(as.numeric(summary(fit_lm_s2)$sigma^2), digits = 8), "\n"
)
cat("dIndependent_Normal_Gamma: Prior_Setup$dispersion =", format(p4_ing$dispersion, digits = 8), "\n")
cat("dNormal_Gamma (n_prior-p): Prior_Setup$dispersion =", format(p4_ng$dispersion, digits = 8), "\n")

cat("\n======== Section 4: posterior vcov(lmb) =========\n")
cat("dIndependent_Normal_Gamma (Post_mean):\n")
print(vcov(fit_p4_ing))
cat("\ndNormal_Gamma (n_prior-p, Post_mean):\n")
print(vcov(fit_p4_ng))

cat("\n======== Section 4: ratios vcov(lmb) / vcov(lm_s2) =========\n")
R4_ing <- vcov(fit_p4_ing) / V_lm_s2
R4_ng <- vcov(fit_p4_ng) / V_lm_s2
cat("dIndependent_Normal_Gamma:\n")
print(R4_ing)
cat("dNormal_Gamma (n_prior-p):\n")
print(R4_ng)
cat(sprintf(
  "Mean diagonal ratio vs lm:  ING %5.3f;  NG (n_prior-p) %5.3f\n",
  mean(diag(R4_ing)),
  mean(diag(R4_ng))
))
R4_ing_s <- vcov(fit_p4_ing) / V_lm_scaled_s2
R4_ng_s <- vcov(fit_p4_ng) / V_lm_scaled_s2
cat("\nSection 4: vcov(lmb) / scaled vcov(lm_s2):\n")
cat(sprintf(
  "Mean diagonal vs scaled:  ING %5.3f;  NG (n_prior-p) %5.3f\n",
  mean(diag(R4_ing_s)),
  mean(diag(R4_ng_s))
))

cat("\nSection 4: congruence vs vcov(lm_s2) / scaled:\n")
cat(congruence_line(vcov(fit_p4_ing), "ING Post_mean", V_lm_s2, "vcov(lm_s2)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p4_ng), "NG Post_mean", V_lm_s2, "vcov(lm_s2)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p4_ing), "ING Post_mean", V_lm_scaled_s2, "scaled vcov(lm_s2)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p4_ng), "NG Post_mean", V_lm_scaled_s2, "scaled vcov(lm_s2)"), "\n", sep = "")

cat("\n")
cat("######################################################################\n")
cat("# Summary table: vcov ratios, posterior dispersion & precision vs lm            #\n")
cat("#   vs_lm  = mean(diag(vcov(lmb) / vcov(lm)))  [Sec 1,3: n=20 lm; 2,4: lm_s2] #\n")
cat("#   n_over_n_plus_nprior = n/(n+n_prior), same n,prior scale as Prior_Setup pwt   #\n")
cat("#   prior_setup_dispersion = Prior_Setup()$dispersion (point estimate for sigma^2)   #\n")
cat("#   vs_scaled = vcov ratio vs (1-pwt)*vcov(lm); note (1-pwt) = n/(n+n_prior) here #\n")
cat("#   target_df_denom = demo goal: n-p (Sec 1,3) vs n (Sec 2,4), not lm$df       #\n")
cat("#   post_mean_disp = mean(lmb$dispersion) - posterior mean sigma^2 (MC)       #\n")
cat("#   disp_vs_lm = post_mean_disp / summary(lm)$sigma^2  (posterior / lm variance) #\n")
cat("#   post_mean_prec = E_post[data precision] = mean(1/dispersion draws), NOT 1/E[disp] #\n")
cat("#   prec_vs_lm = post_mean_prec / (1/dispersion_lm); dispersion_lm = summary(lm)$sigma^2 #\n")
cat("#     Denominator is lm data precision; do NOT use 1/disp_vs_lm.                    #\n")
cat("#     Sec 1,3: fit_lm; Sec 2,4: fit_lm_s2.                                         #\n")
cat("######################################################################\n")
cat(
  "\nPost_mean optimizes the fragment fixed point, not closeness of marginal lmb() vcov\n",
  "to lm; column vs_lm can sit farther from 1 than OLS_mean rows (see header).\n\n",
  sep = ""
)

d_lm_n20 <- as.numeric(summary(fit_lm)$sigma^2)
d_lm_rep <- as.numeric(summary(fit_lm_s2)$sigma^2)
## dispersion_lm = sigma^2_hat from lm (RSS/(n-p)); precision_lm = 1/dispersion_lm
dispersion_lm_tbl <- c(
  d_lm_n20, d_lm_n20, d_lm_rep, d_lm_rep,
  d_lm_n20, d_lm_n20, d_lm_rep, d_lm_rep
)
precision_lm_tbl <- 1 / dispersion_lm_tbl

n_eff_sec1 <- as.integer(stats::nobs(fit_lm))
n_prior_sec1_impl <- (pwt_sec1 / (1 - pwt_sec1)) * n_eff_sec1
n_over_nnp_sec1 <- n_eff_sec1 / (n_eff_sec1 + n_prior_sec1_impl)
n_over_nnp_sec2 <- n_eff_s2 / (n_eff_s2 + n_prior_sec2)
n_over_nnp_tbl <- c(
  n_over_nnp_sec1, n_over_nnp_sec1,
  n_over_nnp_sec2, n_over_nnp_sec2,
  n_over_nnp_sec1, n_over_nnp_sec1,
  n_over_nnp_sec2, n_over_nnp_sec2
)

ratio_summary <- data.frame(
  section = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L),
  disp_type = c(
    "OLS_mean", "OLS_mean", "OLS_mean", "OLS_mean",
    "Post_mean", "Post_mean", "Post_mean", "Post_mean"
  ),
  pfamily = c(
    "NG (shape_df n_prior)",
    "ING (n_prior+p)",
    "ING (shape_df n_prior)",
    "NG (n_prior-p)",
    "NG (shape_df n_prior)",
    "ING (n_prior+p)",
    "ING (shape_df n_prior)",
    "NG (n_prior-p)"
  ),
  target_df_denom = c("n-p", "n-p", "n", "n", "n-p", "n-p", "n", "n"),
  vs_lm = c(
    mean(diag(R1_ng)),
    mean(diag(R1_ing)),
    mean(diag(R2_ing)),
    mean(diag(R2_ng)),
    mean(diag(R3_ng)),
    mean(diag(R3_ing)),
    mean(diag(R4_ing)),
    mean(diag(R4_ng))
  ),
  n_over_n_plus_nprior = n_over_nnp_tbl,
  prior_setup_dispersion = c(
    p1_ng$dispersion,
    p1_ing$dispersion,
    p2_ing$dispersion,
    p2_ng$dispersion,
    p3_ng$dispersion,
    p3_ing$dispersion,
    p4_ing$dispersion,
    p4_ng$dispersion
  ),
  vs_scaled = c(
    mean(diag(R1_ng_s)),
    mean(diag(R1_ing_s)),
    mean(diag(R2_ing_s)),
    mean(diag(R2_ng_s)),
    mean(diag(R3_ng_s)),
    mean(diag(R3_ing_s)),
    mean(diag(R4_ing_s)),
    mean(diag(R4_ng_s))
  ),
  post_mean_disp = c(
    post_mean_lmb_dispersion(fit_p1_ng),
    post_mean_lmb_dispersion(fit_p1_ing),
    post_mean_lmb_dispersion(fit_p2_ing),
    post_mean_lmb_dispersion(fit_p2_ng),
    post_mean_lmb_dispersion(fit_p3_ng),
    post_mean_lmb_dispersion(fit_p3_ing),
    post_mean_lmb_dispersion(fit_p4_ing),
    post_mean_lmb_dispersion(fit_p4_ng)
  ),
  post_mean_prec = c(
    post_mean_lmb_precision(fit_p1_ng),
    post_mean_lmb_precision(fit_p1_ing),
    post_mean_lmb_precision(fit_p2_ing),
    post_mean_lmb_precision(fit_p2_ng),
    post_mean_lmb_precision(fit_p3_ng),
    post_mean_lmb_precision(fit_p3_ing),
    post_mean_lmb_precision(fit_p4_ing),
    post_mean_lmb_precision(fit_p4_ng)
  ),
  stringsAsFactors = FALSE
)
ratio_summary$vs_lm <- round(ratio_summary$vs_lm, 4)
ratio_summary$n_over_n_plus_nprior <- round(ratio_summary$n_over_n_plus_nprior, 6)
ratio_summary$prior_setup_dispersion <- round(ratio_summary$prior_setup_dispersion, 6)
ratio_summary$vs_scaled <- round(ratio_summary$vs_scaled, 4)
ratio_summary$disp_vs_lm <- round(ratio_summary$post_mean_disp / dispersion_lm_tbl, 4)
ratio_summary$prec_vs_lm <- round(ratio_summary$post_mean_prec / precision_lm_tbl, 4)
ratio_summary$post_mean_disp <- round(ratio_summary$post_mean_disp, 6)
ratio_summary$post_mean_prec <- round(ratio_summary$post_mean_prec, 6)
ratio_summary <- ratio_summary[c(
  "section", "disp_type", "pfamily", "target_df_denom",
  "vs_lm", "n_over_n_plus_nprior", "prior_setup_dispersion", "vs_scaled",
  "post_mean_disp", "disp_vs_lm",
  "post_mean_prec", "prec_vs_lm"
)]
print(ratio_summary, row.names = FALSE)
cat(
  "\nReference dispersion_lm = summary(lm)$sigma^2 (disp_vs_lm denominator):  n=20: ",
  format(d_lm_n20, digits = 10),
  " ; replicated: ", format(d_lm_rep, digits = 10), "\n",
  "Reference 1/dispersion_lm = lm data precision (prec_vs_lm denominator):  n=20: ",
  format(precision_lm_tbl[1L], digits = 10),
  " ; replicated: ", format(precision_lm_tbl[3L], digits = 10), "\n",
  sep = ""
)
cat("\nSee ?Prior_Setup, argument shape_df.\n")

invisible(NULL)
