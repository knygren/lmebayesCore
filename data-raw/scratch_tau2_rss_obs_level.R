## RSS from VarCorr: group-level vs observation-level (random slope)
pkgload::load_all(".", quiet = TRUE)
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
gl <- levels(design$groups)
J <- length(gl)
n_obs <- nrow(dat)

lmer_fit <- lmer(form, data = dat, REML = TRUE)
vc <- extract_mer_variance_components(lmer_fit, design$re_coef_names)
lmer_fixef <- fixef(lmer_fit)
lmer_coef <- coef(lmer_fit)[["school_id"]]
lmer_ranef <- ranef(lmer_fit)[["school_id"]]

Z <- design$Z
groups <- design$groups
g_idx <- match(as.character(groups), gl)

cat(sprintf("n_obs = %d, J = %d\n\n", n_obs, J))

for (k in design$re_coef_names) {
  tau2 <- vc$vcov_re[[k]]
  z_k <- Z[, k]
  b_j <- lmer_coef[gl, k]
  u_j <- lmer_ranef[gl, k]
  gk <- if (k == "(Intercept)") {
    lmer_fixef[colnames(design$X_hyper[[k]])]
  } else if (k == "distracted_ppvt") {
    setNames(lmer_fixef["distracted_ppvt"], "(Intercept)")
  } else {
    setNames(c(
      lmer_fixef["distracted_a1"],
      lmer_fixef["free_reduced_lunch:distracted_a1"]
    ), colnames(design$X_hyper[[k]]))
  }
  Xk <- design$X_hyper[[k]]
  eta_j <- as.numeric(Xk[gl, , drop = FALSE] %*% gk)

  ## Group-level RSS (current code)
  rss_grp <- sum((b_j - eta_j)^2)

  ## Observation-level: sum_i (z_ik * u_{g(i)})^2  [BLUP on RE scale]
  u_obs <- u_j[g_idx]
  rss_obs_u <- sum((z_k * u_obs)^2)

  ## Observation-level with total coef b per group: z * (eta + u) - z*eta = z*u
  rss_obs_same <- rss_obs_u

  ## Weighted by z^2 per school: sum_j u_j^2 * sum_{i in j} z_ij^2
  ss_z_by_school <- tapply(z_k^2, groups, sum)[gl]
  rss_weighted_u <- sum(u_j^2 * ss_z_by_school)

  ## If VarCorr = tau2, "equivalent" RSS at obs level (rough): tau2 * sum(z^2)
  rss_from_vc_obs <- tau2 * sum(z_k^2)

  ## Group J*tau2
  rss_j_tau2 <- J * tau2

  cat(sprintf("=== %s  (VarCorr = %.2f) ===\n", k, tau2))
  cat(sprintf("  RSS_grp  = sum_j (b_j - eta_j)^2           = %10.1f\n", rss_grp))
  cat(sprintf("  RSS_obs  = sum_i (z_ik * u_j)^2            = %10.1f\n", rss_obs_u))
  cat(sprintf("  RSS_wt   = sum_j u_j^2 * sum_{i in j} z^2  = %10.1f\n", rss_weighted_u))
  cat(sprintf("  sum(z^2) = %10.1f\n", sum(z_k^2)))
  cat(sprintf("  J*tau2   = %10.1f\n", rss_j_tau2))
  cat(sprintf("  tau2*sum(z^2) (VarCorr-implied obs RSS)    = %10.1f\n", rss_from_vc_obs))
  cat(sprintf("  RSS_grp / (J*tau2) = %.3f\n", rss_grp / rss_j_tau2))
  cat(sprintf("  RSS_obs / (tau2*sum(z^2)) = %.3f\n\n", rss_obs_u / rss_from_vc_obs))
}

## ppvt detail: effective df
k <- "distracted_ppvt"
tau2 <- vc$vcov_re[[k]]
z <- Z[, k]
u_obs <- lmer_ranef[match(groups, gl), k]
cat("ppvt: sum(z^2) =", sum(z^2), ", mean per obs =", mean(z^2), "\n")
cat("ppvt: var(u_obs*z) =", var(u_obs * z), "\n")
cat("pppt: sum((z*u)^2)/sum(z^2) =", sum((z * u_obs)^2) / sum(z^2), " cf VarCorr", tau2, "\n")
