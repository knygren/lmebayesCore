## TEMPORARY / SCRATCH -- investigation only, not package code.
##
## Per-group log-concavity ("RSS ellipsoid") test from
## inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 16.3:
##
##   beta_hat_ols_j = (D_j'D_j)^{-1} D_j'y_j
##   RSS_ols_j      = ||y_j - D_j beta_hat_ols_j||^2
##   threshold_j    = 2*r0_j + RSS_ols_j        (r0_j = calibrated Block-1
##                                                Gamma prior RATE for group j;
##                                                the shape a0_j drops out of
##                                                the final criterion)
##   q_j(beta)      = (beta - beta_hat_ols_j)' D_j'D_j (beta - beta_hat_ols_j)
##
## Locally log-concave (inside the ellipsoid) iff q_j(beta) <= threshold_j.
##
## Reports, per group: whether the POSTERIOR MEAN of beta_j is inside/outside,
## and what percentage of the main-stage MCMC draws of beta_j fall outside.

.tmp_rss_ellipsoid_test <- function(fit, D, y, group, group_name, re_coef_names,
                                     rate_group, group_levels = levels(group)) {
  beta_bar <- lmebayesCore:::.lmebayes_posterior_mean_group_coef(
    fit, group_name, re_coef_names
  )
  co <- fit$coefficients
  group_chr <- as.character(group)

  rows <- lapply(group_levels, function(lev) {
    idx <- which(group_chr == lev)
    D_j <- D[idx, re_coef_names, drop = FALSE]
    y_j <- y[idx]
    DtD <- crossprod(D_j)
    beta_ols <- as.vector(solve(DtD, crossprod(D_j, y_j)))
    RSS_ols  <- sum((y_j - D_j %*% beta_ols)^2)
    thresh   <- 2 * rate_group[[lev]] + RSS_ols

    draws <- as.matrix(co[co[[group_name]] == lev, re_coef_names, drop = FALSE])
    diffs <- sweep(draws, 2, beta_ols, "-")
    q_draws <- rowSums((diffs %*% DtD) * diffs)

    d_mean <- beta_bar[[lev]] - beta_ols
    q_mean <- as.numeric(t(d_mean) %*% DtD %*% d_mean)

    data.frame(
      group = lev,
      n_j = length(idx),
      RSS_ols = RSS_ols,
      threshold = thresh,
      q_mean = q_mean,
      inside_at_mean = q_mean <= thresh,
      pct_draws_outside = 100 * mean(q_draws > thresh)
    )
  })
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## Usage: with an existing fit already in your session (Ex_13/Ex_13b/Ex_14's
## variable names: fit, design, grp, re_names, rate_group), just call:
##
##   tab <- .tmp_rss_ellipsoid_test(
##     fit = fit, D = design$D, y = design$y, group = grp,
##     group_name = design$group_name, re_coef_names = re_names,
##     rate_group = rate_group
##   )
##   print(tab, row.names = FALSE, digits = 4)
##   cat(sprintf("Groups with posterior mean OUTSIDE the ellipsoid: %d / %d\n",
##               sum(!tab$inside_at_mean), nrow(tab)))
##   cat(sprintf("Mean pct of draws outside (across groups): %.2f%%\n",
##               mean(tab$pct_draws_outside)))
##
## No demo re-run needed -- this file only *defines* .tmp_rss_ellipsoid_test()
## above; sourcing it does not execute anything else.
## ---------------------------------------------------------------------------
