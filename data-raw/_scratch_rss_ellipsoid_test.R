## TEMPORARY / SCRATCH -- investigation only, not package code.
##
## Per-group log-concavity ("RSS ellipsoid") test from
## inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 16.3:
##
##   beta_hat_ols_j = (D_j'D_j)^{-1} D_j'y_j
##   RSS_ols_j      = ||y_j - D_j beta_hat_ols_j||^2
##   threshold_j    = 2*r0_j + RSS_ols_j        (r0_j = calibrated Block-1
##                                                Gamma prior RATE for group j)
##   q_j(beta)      = (beta - beta_hat_ols_j)' D_j'D_j (beta - beta_hat_ols_j)
##
## Locally log-concave (inside the ellipsoid) iff q_j(beta) <= threshold_j.
##
## Reports, per group: whether the POSTERIOR MEAN of beta_j is inside/outside,
## what percentage of the main-stage MCMC draws of beta_j fall outside, and
## (new) an observed-vs-expected comparison against a no-shrinkage Bayesian-t
## baseline, per inst/multivariate-t-log-concavity.md's "Confidence-level
## interpretation via the Scheffe region" -- generalized from the classical
## n-p residual df to nu_j = 2*a0_j + n_j - p_re (a0_j = Block-1 ING Gamma
## prior SHAPE for group j, folding in the Gamma prior's own pseudo-df):
##
##   alpha_j  = P(F_{p_re,nu_j} > nu_j/p_re)   -- tail fraction of beta_j's
##                own (unshrunk) t-posterior that the ellipsoid boundary
##                excludes, i.e. the "expected pct outside" if group j had
##                no pull from Lambda/W_j*gamma at all.
##   obs_n_j  = count of main-stage draws with q_j(beta_j) > threshold_j
##   exp_n_j  = alpha_j * n_draws_j
##   ratio_j  = obs_n_j / exp_n_j              (descriptive only -- unstable
##                near alpha_j ~ 0, see z_j/p_value_j below)
##   z_j      = (obs_n_j - exp_n_j) / sqrt(n_draws_j*alpha_j*(1-alpha_j))
##   p_value_j = P(Binomial(n_draws_j, alpha_j) >= obs_n_j)   (one-sided;
##                small p_value_j flags a group with far more draws outside
##                the ellipsoid than its own data+prior alone would predict)
##
## z_j/p_value_j treat the n_draws_j "outside" indicators as iid, which MCMC
## autocorrelation violates -- fine for a first-pass ranking/screening tool,
## not a calibrated p-value.

.tmp_rss_ellipsoid_test <- function(fit, D, y, group, group_name, re_coef_names,
                                     shape_group, rate_group,
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

    draws <- as.matrix(co[co[[group_name]] == lev, re_coef_names, drop = FALSE])
    n_draws_j <- nrow(draws)
    diffs <- sweep(draws, 2, beta_ols, "-")
    q_draws <- rowSums((diffs %*% DtD) * diffs)
    obs_n <- sum(q_draws > thresh)

    d_mean <- beta_bar[[lev]] - beta_ols
    q_mean <- as.numeric(t(d_mean) %*% DtD %*% d_mean)

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
## variable names: fit, design, grp, re_names, shape_group, rate_group), just
## call:
##
##   tab <- .tmp_rss_ellipsoid_test(
##     fit = fit, D = design$D, y = design$y, group = grp,
##     group_name = design$group_name, re_coef_names = re_names,
##     shape_group = shape_group, rate_group = rate_group
##   )
##   print(tab[order(tab$p_value), ], row.names = FALSE, digits = 4)
##   cat(sprintf("Groups with posterior mean OUTSIDE the ellipsoid: %d / %d\n",
##               sum(!tab$inside_at_mean), nrow(tab)))
##   cat(sprintf("Mean pct of draws outside (across groups): %.2f%%\n",
##               mean(tab$pct_draws_outside)))
##
## No demo re-run needed -- this file only *defines* .tmp_rss_ellipsoid_test()
## above; sourcing it does not execute anything else.
## ---------------------------------------------------------------------------
