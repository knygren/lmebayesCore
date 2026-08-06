## TEMPORARY / SCRATCH -- investigation only.
## Checks whether "0.000% outside for every group" (both pre-run calibration
## prediction and post-run marginal ellipsoid check on Ex_13b) reflects a
## genuine (if surprising) consequence of max_disp_perc_measurement = 0.8
## shrinking Var_t[Omega_j] enough that K_j = E_t[Omega_j]/Var_t[Omega_j] is
## just very large everywhere, or a bug making K_j blow up / q_j collapse.

devtools::load_all(".", quiet = TRUE)

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

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

design_all <- model_setup(form_lmer, data = dat)
full_rank_schools <- names(design_all$re_rank)[design_all$re_rank]
dat <- subset(dat, school_id %in% full_rank_schools)
dat$school_id <- droplevels(dat$school_id)

ps <- Prior_Setup_GLMM(
  form_lmer,
  data            = dat,
  pwt             = 0.01,
  dispformula     = ~school_id,
  max_disp_perc_measurement = 0.8,
  pwt_measurement = 0.1,
  alpha_target_measurement  = 0.01
)

cat("\n=== pwt_measurement_calibration (full precision) ===\n\n")
calib <- ps$pwt_measurement_calibration
print(calib, row.names = FALSE, digits = 10)

cat(sprintf(
  "\nGroups where pct_outside_before <= alpha_target (0.01) already at floor w=0.1: %d / %d\n",
  sum(calib$pct_outside_before <= 1), nrow(calib)
))
cat(sprintf(
  "Range of pct_outside_before (%%%%): [%.6g, %.6g]\n",
  min(calib$pct_outside_before), max(calib$pct_outside_before)
))
cat(sprintf(
  "Range of pct_outside_after  (%%%%): [%.6g, %.6g]\n",
  min(calib$pct_outside_after), max(calib$pct_outside_after)
))
cat(sprintf(
  "clipped_at_ceiling: %d groups\n", sum(calib$clipped_at_ceiling)
))
cat(sprintf(
  "floor_binds: %d groups\n", sum(calib$floor_binds)
))

## ---------------------------------------------------------------------------
## Spot-check ONE group's actual K_j vs q_j distribution directly (bypassing
## the calibration search entirely) to see whether K_j is merely "large" or
## pathologically huge/Inf.
## ---------------------------------------------------------------------------
group_levels <- levels(dat$school_id)
lev <- group_levels[[1]]
re_names <- design_all$re_coef_names
idx <- which(as.character(dat$school_id) == lev)

design <- model_setup(form_lmer, data = dat)
D_j <- design$D[which(as.character(design$group) == lev), re_names, drop = FALSE]
y_j <- design$y[which(as.character(design$group) == lev)]
n_j <- length(y_j)
DtD <- crossprod(D_j)
beta_ols <- as.vector(solve(DtD, crossprod(D_j, y_j)))
RSS_ols  <- sum((y_j - D_j %*% beta_ols)^2)

ing_j <- ps$ing_prior_measurement_group[[lev]]
a0 <- ing_j$shape_ING
r0 <- ing_j$rate
disp_lower <- ing_j$disp_lower
disp_upper <- ing_j$disp_upper
omega_L <- 1 / disp_upper
omega_U <- 1 / disp_lower
shape_j <- a0 + n_j / 2

cat(sprintf(
  "\n=== Spot-check group '%s': n_j=%d, a0=%.4f, r0=%.4f, shape_j=%.4f ===\n",
  lev, n_j, a0, r0, shape_j
))
cat(sprintf("disp_lower=%.6g, disp_upper=%.6g, omega_L=%.6g, omega_U=%.6g\n",
            disp_lower, disp_upper, omega_L, omega_U))

## Simulate beta_j draws from its own Gaussian full conditional at a plug-in
## Omega_hat_j to get a realistic q_j spread, then evaluate K_j at each.
Sigma_ranef <- ps$Sigma_ranef
Psi_inv <- solve(Sigma_ranef)
Omega_hat_j <- r0 / (a0 - 1) # prior mean precision proxy: 1/sigma2_hat roughly
sigma2_hat_j <- ing_j$sigma2_hat
Omega_hat_j <- 1 / sigma2_hat_j

Sigma_j_inv <- Omega_hat_j * DtD + Psi_inv
Sigma_j <- solve(Sigma_j_inv)
set.seed(1)
n_sim <- 20000L
p_re <- length(re_names)
L <- chol(Sigma_j)
Z <- matrix(rnorm(n_sim * p_re), nrow = n_sim)
draws <- sweep(Z %*% L, 2, beta_ols, "+")  # centered loosely at beta_ols for spot-check
diffs <- sweep(draws, 2, beta_ols, "-")
q_draws <- rowSums((diffs %*% DtD) * diffs)

rate_draws <- r0 + 0.5 * (q_draws + RSS_ols)
mom <- lmebayesCore:::.two_block_truncated_omega_moments(shape_j, rate_draws, omega_L, omega_U)
K_draws <- mom$E_omega / mom$Var_omega

cat(sprintf("q_draws:  mean=%.4g, sd=%.4g, range=[%.4g, %.4g]\n",
            mean(q_draws), sd(q_draws), min(q_draws), max(q_draws)))
cat(sprintf("K_draws:  mean=%.4g, sd=%.4g, range=[%.4g, %.4g]\n",
            mean(K_draws), sd(K_draws), min(K_draws), max(K_draws)))
cat(sprintf("E_omega:  mean=%.6g, range=[%.6g, %.6g]\n",
            mean(mom$E_omega), min(mom$E_omega), max(mom$E_omega)))
cat(sprintf("Var_omega: mean=%.6g, range=[%.6g, %.6g]\n",
            mean(mom$Var_omega), min(mom$Var_omega), max(mom$Var_omega)))
cat(sprintf("pct q_draws > K_draws: %.4f%%\n", 100 * mean(q_draws > K_draws)))

## Untruncated Var(Omega) for comparison (Gamma(shape_j, rate_draws) without
## truncation): Var = shape/rate^2. Ratio truncated/untruncated shows how
## much the window is shrinking the variance.
Var_untrunc <- shape_j / rate_draws^2
cat(sprintf(
  "Var_omega (truncated) / Var (untruncated) ratio: mean=%.6g, range=[%.6g, %.6g]\n",
  mean(mom$Var_omega / Var_untrunc), min(mom$Var_omega / Var_untrunc), max(mom$Var_omega / Var_untrunc)
))
