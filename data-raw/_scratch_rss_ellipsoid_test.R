## TEMPORARY / SCRATCH -- investigation only, not package code.
##
## Per-group log-concavity ("RSS ellipsoid") test from
## inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 16.3/16.6:
##
##   beta_hat_ols_j = (D_j'D_j)^{-1} D_j'y_j
##   RSS_ols_j      = ||y_j - D_j beta_hat_ols_j||^2
##   q_j(beta)      = (beta - beta_hat_ols_j)' D_j'D_j (beta - beta_hat_ols_j)
##                    = e_j(beta)'e_j(beta) - RSS_ols_j  (Pythagorean OLS split)
##
## Two versions of the boundary are reported:
##
##  UNTRUNCATED (Section 16.3), a FIXED per-group number, kept as `threshold`/
##  `pct_draws_outside`/etc. for continuity with earlier runs:
##
##    threshold_j = 2*r0_j + RSS_ols_j   (r0_j = calibrated Block-1 Gamma
##                                         prior RATE for group j)
##
##  EXACT / truncation-aware (Section 16.6), a PER-DRAW number, since the
##  actual sampled Omega_j is Gamma(shape_j, rate_j(beta)) truncated to
##  [omega_L_j, omega_U_j] = [1/disp_upper_j, 1/disp_lower_j], and the
##  truncated moments E_t[Omega_j]/Var_t[Omega_j] depend on that draw's own
##  rate_j(beta) = r0_j + 0.5*e_j(beta)'e_j(beta):
##
##    threshold_exact_j(beta) = E_t[Omega_j] / Var_t[Omega_j]
##
##  (via lmebayesCore:::.two_block_truncated_omega_moments(), vectorized over
##  draws for a fixed group). Locally log-concave (inside the exact-Hessian
##  ellipsoid) iff q_j(beta) <= threshold_exact_j(beta); this is stricter
##  book-keeping than the untruncated version but algebraically identical to
##  the "raw H_j alone" criterion used by
##  .tmp_lambda_star_marginal_over_draws()'s h_violates_count_exact.
##
## Reports, per group: whether the POSTERIOR MEAN of beta_j is inside/outside
## (both versions), what percentage of the main-stage MCMC draws of beta_j
## fall outside (both versions), and an observed-vs-expected comparison
## against a no-shrinkage Bayesian-t baseline, per
## inst/multivariate-t-log-concavity.md's "Confidence-level interpretation
## via the Scheffe region" -- generalized from the classical n-p residual df
## to nu_j = 2*a0_j + n_j - p_re (a0_j = Block-1 ING Gamma prior SHAPE for
## group j, folding in the Gamma prior's own pseudo-df):
##
##   alpha_j  = P(F_{p_re,nu_j} > nu_j/p_re)   -- tail fraction of beta_j's
##                own (unshrunk) t-posterior that the UNTRUNCATED ellipsoid
##                boundary excludes, i.e. the "expected pct outside" if
##                group j had no pull from Lambda/W_j*gamma at all. This
##                F-pivot baseline is untruncated-only (Section 16.3's a0_j
##                cancels out of threshold_j but NOT out of threshold_exact_j
##                -- see Section 16.6/9 for why the exact criterion has no
##                comparably simple closed-form baseline); it is still shown
##                alongside the exact obs_n/pct columns as a conservative
##                reference, not a matched null.
##   obs_n_j  = count of main-stage draws with q_j(beta_j) > threshold_j
##                (UNTRUNCATED)
##   exp_n_j  = alpha_j * n_draws_j
##   ratio_j  = obs_n_j / exp_n_j              (descriptive only -- unstable
##                near alpha_j ~ 0, see z_j/p_value_j below)
##   z_j      = (obs_n_j - exp_n_j) / sqrt(n_draws_j*alpha_j*(1-alpha_j))
##   p_value_j = P(Binomial(n_draws_j, alpha_j) >= obs_n_j)   (one-sided;
##                small p_value_j flags a group with far more draws outside
##                the UNTRUNCATED ellipsoid than its own data+prior alone
##                would predict)
##
## z_j/p_value_j treat the n_draws_j "outside" indicators as iid, which MCMC
## autocorrelation violates -- fine for a first-pass ranking/screening tool,
## not a calibrated p-value.

## =============================================================================
## .tmp_rss_ellipsoid_test_marginal(): THE ellipsoid check (Omega integrated
## out exactly, Section 16.6) -- add this new section BEFORE the older,
## untruncated/no-shrinkage-baseline test below.
##
## The whole point of calibrating disp_lower_j/disp_upper_j
## (Prior_Setup_lmebayes()'s alpha_target_measurement search,
## R/pwt_measurement_group_calibration.R) is to keep the fraction of
## log-concavity violations below a target. That search is a PRE-RUN
## approximation: it simulates q_j from a plug-in NONCENTRAL GAUSSIAN full
## conditional (Omega_j, gamma, Lambda all fixed at reference point
## estimates -- see .lmebayes_calibrate_pwt_measurement_group()). This
## function is the POST-RUN check on whether that prediction actually held
## up against the real posterior: for each of the n main-stage draws
## (independent parallel chains, NOT a single autocorrelated chain -- see
## m_convergence/Theorem-3 calibration), it uses ONLY that draw's own beta_j
## (via q_j) to test
##
##   q_j(beta_j) <= K_j(beta_j) := E_t[Omega_j | beta_j] / Var_t[Omega_j | beta_j]
##
## i.e. Omega_j analytically integrated out of the marginal Hessian
## H_j(beta_j) = E_t[Omega_j]*D_j'D_j - Var_t[Omega_j]*(D_j'e_j)(D_j'e_j)'
## (Section 16.6) -- NOT the realized dispersion_ranef draw for that sweep
## (conditioning on a specific, realized Omega_j makes H_j trivially PD for
## every beta_j, which is why Omega_j must be integrated out analytically
## instead: see the "why not just use the sampled Omega_j" discussion this
## section grew out of). A draw is flagged "degenerate" (excluded from both
## the numerator and denominator of pct_draws_outside below) only when the
## truncated-moment calculation itself is numerically unstable for that
## draw's own q_j -- see .two_block_truncated_omega_moments()'s dP0_floor
## guard -- not because any Omega_j draw ever left [omega_L, omega_U] (the
## sampler enforces that truncation directly; it never happens).
##
## Pass `pwt_measurement_calibration = ps$pwt_measurement_calibration` (the
## data.frame Prior_Setup_lmebayes() stores whenever alpha_target_measurement
## calibration ran) to attach each group's PRE-RUN predicted
## `pct_outside_after` alongside the POST-RUN empirical `pct_draws_outside`
## computed here, for direct comparison.
## =============================================================================
.tmp_rss_ellipsoid_test_marginal <- function(fit, D, y, group, group_name, re_coef_names,
                                              shape_group, rate_group,
                                              disp_lower_group, disp_upper_group,
                                              group_levels = levels(group),
                                              pwt_measurement_calibration = NULL) {
  co <- fit$coefficients
  group_chr <- as.character(group)

  rows <- lapply(group_levels, function(lev) {
    idx <- which(group_chr == lev)
    D_j <- D[idx, re_coef_names, drop = FALSE]
    y_j <- y[idx]
    n_j <- length(idx)
    DtD <- crossprod(D_j)
    beta_ols <- as.vector(solve(DtD, crossprod(D_j, y_j)))
    RSS_ols  <- sum((y_j - D_j %*% beta_ols)^2)

    a0 <- shape_group[[lev]]
    r0 <- rate_group[[lev]]
    omega_L  <- 1 / disp_upper_group[[lev]]
    omega_U  <- 1 / disp_lower_group[[lev]]
    shape_j  <- a0 + n_j / 2

    draws <- as.matrix(co[co[[group_name]] == lev, re_coef_names, drop = FALSE])
    n_draws_j <- nrow(draws)
    diffs <- sweep(draws, 2, beta_ols, "-")
    q_draws <- rowSums((diffs %*% DtD) * diffs)

    ## rate_j(beta) = r0 + 0.5*e_j(beta)'e_j(beta) = r0 + 0.5*(q_j(beta) + RSS_ols_j)
    rate_draws <- r0 + 0.5 * (q_draws + RSS_ols)
    mom_draws  <- lmebayesCore:::.two_block_truncated_omega_moments(
      shape_j, rate_draws, omega_L, omega_U
    )
    K_draws      <- mom_draws$E_omega / mom_draws$Var_omega
    n_degenerate <- sum(!mom_draws$ok)
    n_eval       <- sum(mom_draws$ok)
    obs_n        <- sum(mom_draws$ok & (q_draws > K_draws))

    data.frame(
      group             = lev,
      n_j               = n_j,
      n_draws           = n_draws_j,
      n_degenerate      = n_degenerate,
      n_eval            = n_eval,
      obs_n             = obs_n,
      pct_draws_outside = 100 * obs_n / max(n_eval, 1L),
      stringsAsFactors  = FALSE
    )
  })
  tab <- do.call(rbind, rows)

  if (!is.null(pwt_measurement_calibration)) {
    pre <- pwt_measurement_calibration[
      match(tab$group, pwt_measurement_calibration$group),
      "pct_outside_after"
    ]
    tab$pct_outside_prerun <- pre
    tab$diff_vs_prerun      <- tab$pct_draws_outside - tab$pct_outside_prerun
  }
  tab
}

## Usage:
##   tab_marg <- .tmp_rss_ellipsoid_test_marginal(
##     fit = fit, D = design$D, y = design$y, group = grp,
##     group_name = design$group_name, re_coef_names = re_names,
##     shape_group = shape_group, rate_group = rate_group,
##     disp_lower_group = disp_lower_group, disp_upper_group = disp_upper_group,
##     pwt_measurement_calibration = ps$pwt_measurement_calibration
##   )
##   print(tab_marg[order(-tab_marg$pct_draws_outside), ], row.names = FALSE, digits = 4)
##   cat(sprintf(
##     "Mean pct outside -- post-run (real posterior) vs pre-run (calibration prediction): %.3f%% vs %.3f%%\n",
##     mean(tab_marg$pct_draws_outside), mean(tab_marg$pct_outside_prerun, na.rm = TRUE)
##   ))

.tmp_rss_ellipsoid_test <- function(fit, D, y, group, group_name, re_coef_names,
                                     shape_group, rate_group,
                                     disp_lower_group, disp_upper_group,
                                     group_levels = levels(group)) {
  beta_bar <- lmebayesCore:::.lmebayes_posterior_mean_group_coef(
    fit, group_name, re_coef_names
  )
  co <- fit$coefficients
  group_chr <- as.character(group)
  p_re <- length(re_coef_names)

  rows <- lapply(group_levels, function(lev) {
    idx <- which(group_chr == lev)
    D_j <- D[idx, re_coef_names, drop = FALSE]
    y_j <- y[idx]
    DtD <- crossprod(D_j)
    beta_ols <- as.vector(solve(DtD, crossprod(D_j, y_j)))
    RSS_ols  <- sum((y_j - D_j %*% beta_ols)^2)
    a0 <- shape_group[[lev]]
    r0 <- rate_group[[lev]]
    thresh   <- 2 * r0 + RSS_ols
    omega_L  <- 1 / disp_upper_group[[lev]]
    omega_U  <- 1 / disp_lower_group[[lev]]
    shape_j  <- a0 + length(idx) / 2

    draws <- as.matrix(co[co[[group_name]] == lev, re_coef_names, drop = FALSE])
    n_draws_j <- nrow(draws)
    diffs <- sweep(draws, 2, beta_ols, "-")
    q_draws <- rowSums((diffs %*% DtD) * diffs)
    obs_n <- sum(q_draws > thresh)

    ## EXACT/truncated: rate_j(beta) = r0 + 0.5*e_j'e_j = r0 + 0.5*(q_j + RSS_ols_j),
    ## vectorized over all n_draws_j draws for this group.
    rate_draws <- r0 + 0.5 * (q_draws + RSS_ols)
    mom_draws  <- lmebayesCore:::.two_block_truncated_omega_moments(
      shape_j, rate_draws, omega_L, omega_U
    )
    thresh_exact_draws <- mom_draws$E_omega / mom_draws$Var_omega
    n_degenerate <- sum(!mom_draws$ok)
    obs_n_exact  <- sum(mom_draws$ok & (q_draws > thresh_exact_draws))
    n_eval_exact <- sum(mom_draws$ok)

    d_mean <- beta_bar[[lev]] - beta_ols
    q_mean <- as.numeric(t(d_mean) %*% DtD %*% d_mean)
    rate_mean <- r0 + 0.5 * (q_mean + RSS_ols)
    mom_mean  <- lmebayesCore:::.two_block_truncated_omega_moments(
      shape_j, rate_mean, omega_L, omega_U
    )
    thresh_exact_mean <- mom_mean$E_omega / mom_mean$Var_omega

    n_j     <- length(idx)
    nu_j    <- 2 * a0 + n_j - p_re
    alpha_j <- stats::pf(nu_j / p_re, p_re, nu_j, lower.tail = FALSE)
    exp_n   <- alpha_j * n_draws_j
    z_j     <- (obs_n - exp_n) / sqrt(n_draws_j * alpha_j * (1 - alpha_j))
    p_j     <- stats::pbinom(obs_n - 1, n_draws_j, alpha_j, lower.tail = FALSE)

    data.frame(
      group = lev,
      n_j = n_j,
      RSS_ols = RSS_ols,
      threshold = thresh,
      q_mean = q_mean,
      inside_at_mean = q_mean <= thresh,
      pct_draws_outside = 100 * obs_n / n_draws_j,
      threshold_exact_mean = thresh_exact_mean,
      inside_at_mean_exact = isTRUE(mom_mean$ok) && q_mean <= thresh_exact_mean,
      pct_draws_outside_exact = 100 * obs_n_exact / max(n_eval_exact, 1L),
      n_degenerate = n_degenerate,
      nu_j = nu_j,
      expected_pct = 100 * alpha_j,
      obs_n = obs_n,
      exp_n = exp_n,
      ratio = ifelse(exp_n > 0, obs_n / exp_n, NA_real_),
      z = z_j,
      p_value = p_j
    )
  })
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## Usage: with an existing fit already in your session (Ex_13/Ex_13b/Ex_13c's
## variable names: fit, design, grp, re_names, shape_group, rate_group,
## disp_lower_group, disp_upper_group), just call:
##
##   tab <- .tmp_rss_ellipsoid_test(
##     fit = fit, D = design$D, y = design$y, group = grp,
##     group_name = design$group_name, re_coef_names = re_names,
##     shape_group = shape_group, rate_group = rate_group,
##     disp_lower_group = disp_lower_group, disp_upper_group = disp_upper_group
##   )
##   print(tab[order(tab$p_value), ], row.names = FALSE, digits = 4)
##   cat(sprintf("Groups with posterior mean OUTSIDE the ellipsoid (untruncated): %d / %d\n",
##               sum(!tab$inside_at_mean), nrow(tab)))
##   cat(sprintf("Groups with posterior mean OUTSIDE the ellipsoid (exact): %d / %d\n",
##               sum(!tab$inside_at_mean_exact), nrow(tab)))
##   cat(sprintf("Mean pct of draws outside, untruncated vs. exact (across groups): %.2f%% vs %.2f%%\n",
##               mean(tab$pct_draws_outside), mean(tab$pct_draws_outside_exact)))
##
## No demo re-run needed -- this file only *defines* .tmp_rss_ellipsoid_test()
## above; sourcing it does not execute anything else.
## ---------------------------------------------------------------------------
