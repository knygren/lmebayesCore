## TEMPORARY / SCRATCH -- investigation only, not package code.
##
## Per-group pwt_measurement large enough to make the pre-run (closed-form
## Gaussian plug-in, i.e. NONCENTRAL -- see inst/omega-ing-marginal-
## multivariate-t.md Section 9) estimate of pct_draws_outside fall at/below a
## target alpha_target, for every group -- as opposed to the earlier (WRONG
## for outlier groups) closed-form recipe that used the central-F "no-
## shrinkage baseline" alpha_j from Section 7-8, which is blind to each
## group's own noncentrality (how far the hierarchy pulls beta_j away from
## beta_hat_ols_j).
##
## Key simplification (Section 9.3, still valid): sigma2_hat_j, hence
## Omega_hat_j, B_j, beta_bar_j, and the whole DISTRIBUTION of
## q_j = (beta_j-beta_hat_ols_j)' D_j'D_j(beta_j-beta_hat_ols_j) is invariant
## to pwt_measurement -- only the BOUNDARY moves as pwt_measurement (hence
## n_prior_j/a0_j/r0_j) changes. So ONE Monte Carlo simulation of q_j's
## distribution per group is enough to reuse across every candidate w_j.
##
## UPDATE (Section 16.6 exact criterion): the boundary is now the EXACT,
## truncation-aware one, q_j > E_t[Omega_j]/Var_t[Omega_j] (a PER-DRAW
## number, since the truncated-Gamma moments depend on that draw's own
## rate_j(beta) = r0_j + 0.5*e_j'e_j), rather than the untruncated FIXED
## threshold 2*r0_j+RSS_ols_j. Under this exact criterion a0_j no longer
## cancels out of the boundary (Section 16.6/9's whole point), so the old
## closed-form single-quantile-lookup inversion for w_j no longer applies;
## it is kept below ONLY as `w_star_noncentral_untrunc`, a legacy/contrast
## column. The new `w_star_exact` column instead 1-D `uniroot()`s over w
## directly against the exact per-draw criterion, still reusing the SAME
## one-time q_draws simulation (only the boundary function of w is
## recomputed per uniroot evaluation -- no re-simulation).
##
## The exact boundary also needs the truncation window
## [omega_L(w), omega_U(w)] = [1/disp_upper(w), 1/disp_lower(w)] at each
## candidate w, which -- per R/mixed_rmerb_helpers.R's
## .lmebayes_calibrate_ing_prior_measurement_group()/
## .lmebayes_ing_prior_quantile_window() -- are just the max_disp_perc/
## (1-max_disp_perc) quantiles of Gamma(a0(w), r0(w)) itself (cheap, no glm
## refit needed); this mirrors the SAME "sigma2_hat_j fixed, only a0_j/r0_j
## move together with n_prior_j(w)" assumption the old closed-form already
## used for the boundary itself, now extended to the truncation window too.
##
## THIS ONLY COMPUTES A VECTOR TO INSPECT -- nothing here feeds back into
## Prior_Setup_lmebayes()/dGamma_list()/the sampler.

## Map hyper-predictor 'nm' of RE component 'k' to fixef(fit_ref)'s name:
## a main effect for intercept-associated predictors, but the OBSERVATION-
## LEVEL INTERACTION TERM for a non-intercept RE component's own hyper-
## predictors -- same mapping Ex_13's own Section 7 table uses.
.tmp_fe_name <- function(k, nm, fe_names) {
  if (identical(k, "(Intercept)") && identical(nm, "(Intercept)")) {
    return("(Intercept)")
  }
  if (identical(nm, "(Intercept)")) return(k)
  if (identical(k, "(Intercept)")) return(nm)
  cand <- c(paste0(nm, ":", k), paste0(k, ":", nm))
  hit <- cand[cand %in% fe_names]
  if (length(hit)) hit[1L] else NA_character_
}

.tmp_group_pwt_measurement_noncentral <- function(ps, design, group, group_levels,
                                                    re_coef_names,
                                                    alpha_target = 0.01,
                                                    w_floor = 0.1,
                                                    n_sim = 200000L,
                                                    seed = 1L) {
  p_re <- length(re_coef_names)

  gamma_hat_raw <- lmebayesCore:::.lmebayes_reference_fixef(ps$fit_ref)
  Sigma_ranef <- as.matrix(ps$Sigma_ranef)
  Psi_inv <- diag(1 / diag(Sigma_ranef))
  dimnames(Psi_inv) <- dimnames(Sigma_ranef)

  gamma_hat_by_component <- stats::setNames(
    lapply(re_coef_names, function(k) {
      Wk <- design$W[[k]]
      nm <- colnames(Wk)
      fe_nm <- vapply(nm, .tmp_fe_name, character(1), k = k, fe_names = names(gamma_hat_raw))
      if (anyNA(fe_nm)) {
        stop("Could not map hyper-predictor(s) for '", k, "' to fixef(fit_ref).", call. = FALSE)
      }
      stats::setNames(unname(gamma_hat_raw[fe_nm]), nm)
    }),
    re_coef_names
  )

  ## Naive CENTRAL-only w_j* (Section 7-8's no-shrinkage baseline) -- kept
  ## only as a contrast column showing what ignoring Delta_j would give.
  g_central <- function(nu) stats::pf(nu / p_re, p_re, nu, lower.tail = FALSE) - alpha_target
  nu_star_central <- stats::uniroot(g_central, lower = p_re + 1e-6, upper = 1e6)$root

  set.seed(seed)
  group_chr <- as.character(group)

  rows <- lapply(group_levels, function(lev) {
    idx <- which(group_chr == lev)
    D_j <- design$D[idx, re_coef_names, drop = FALSE]
    y_j <- design$y[idx]
    n_j <- length(idx)
    DtD <- crossprod(D_j)
    beta_ols <- as.vector(solve(DtD, crossprod(D_j, y_j)))
    RSS_ols  <- sum((y_j - D_j %*% beta_ols)^2)

    ing_j         <- ps$group.ing_prior[[lev]]
    sigma2_hat    <- ing_j$sigma2_hat
    Omega_hat_j   <- 1 / sigma2_hat
    rate_now      <- ing_j$rate           ## r0_j at ps's CURRENT pwt_measurement
    thresh_now    <- 2 * rate_now + RSS_ols
    max_disp_perc <- ing_j$max_disp_perc

    mu_j <- vapply(re_coef_names, function(k) {
      Wk <- design$W[[k]]
      as.numeric(Wk[lev, , drop = FALSE] %*% gamma_hat_by_component[[k]])
    }, numeric(1))

    Sigma_j_inv <- Omega_hat_j * DtD + Psi_inv
    Sigma_j     <- solve(Sigma_j_inv)
    beta_bar_j  <- as.vector(Sigma_j %*% (Omega_hat_j * crossprod(D_j, y_j) + Psi_inv %*% mu_j))

    ## ONE simulation of q_j's distribution -- invariant to pwt_measurement.
    L <- chol(Sigma_j)
    Z <- matrix(rnorm(n_sim * p_re), nrow = n_sim)
    draws <- sweep(Z %*% L, 2, beta_bar_j, "+")
    diffs <- sweep(draws, 2, beta_ols, "-")
    q_draws <- rowSums((diffs %*% DtD) * diffs)
    ete_draws <- q_draws + RSS_ols

    pct_outside_now <- 100 * mean(q_draws > thresh_now)

    ## EXACT criterion at the CURRENT (ps's) window, for comparison with the
    ## untruncated pct_outside_now above.
    omega_L_now <- 1 / ing_j$disp_upper
    omega_U_now <- 1 / ing_j$disp_lower
    shape_j_now <- ing_j$shape_ING + n_j / 2
    rate_j_now  <- rate_now + 0.5 * ete_draws
    mom_now <- lmebayesCore:::.two_block_truncated_omega_moments(
      shape_j_now, rate_j_now, omega_L_now, omega_U_now
    )
    n_ok_now <- sum(mom_now$ok)
    pct_outside_now_exact <- if (n_ok_now > 0L) {
      100 * sum(mom_now$ok & (q_draws > mom_now$E_omega / mom_now$Var_omega)) / n_ok_now
    } else {
      NA_real_
    }

    ## LEGACY closed-form (Section 7-8/9.3, UNTRUNCATED-boundary): quantile-
    ## match to alpha_target, then invert r0 -> a0 -> n_prior -> w. Kept only
    ## as a contrast column; a0_j cancels out of the untruncated boundary,
    ## which is what makes this single-shot inversion possible there.
    q_star   <- as.numeric(stats::quantile(q_draws, probs = 1 - alpha_target, type = 7))
    r0_star  <- (q_star - RSS_ols) / 2
    a0_star  <- r0_star / sigma2_hat + 1
    n_prior_star <- 2 * a0_star - p_re - 1
    w_raw    <- n_prior_star / (n_prior_star + n_j)
    w_star   <- max(0, min(0.5, w_raw))
    clipped  <- w_raw > 0.5
    w_star_floored <- max(w_star, w_floor)

    w_central <- max(0, (nu_star_central - n_j - 1) / (nu_star_central - 1))

    ## EXACT criterion: 1-D uniroot() over w, reusing q_draws/ete_draws.
    ## a0(w)/r0(w) keep sigma2_hat FIXED (same assumption the legacy closed
    ## form already relies on for the boundary) while n_prior(w) = w/(1-w)*n_j
    ## moves a0/r0 together; disp_lower(w)/disp_upper(w) are recomputed from
    ## the SAME a0(w)/r0(w) via .lmebayes_ing_prior_quantile_window(), so the
    ## truncation window tightens/widens consistently with the prior's own
    ## sharpness as w changes.
    pct_outside_exact_for_w <- function(w) {
      n_prior_w <- w / (1 - w) * n_j
      a0_w <- (n_prior_w + p_re + 1) / 2
      r0_w <- sigma2_hat * (a0_w - 1)
      win_w <- lmebayesCore:::.lmebayes_ing_prior_quantile_window(a0_w, r0_w, max_disp_perc)
      omega_L_w <- 1 / win_w$disp_upper
      omega_U_w <- 1 / win_w$disp_lower
      shape_j_w <- a0_w + n_j / 2
      rate_j_w  <- r0_w + 0.5 * ete_draws
      mom_w <- lmebayesCore:::.two_block_truncated_omega_moments(
        shape_j_w, rate_j_w, omega_L_w, omega_U_w
      )
      n_ok <- sum(mom_w$ok)
      if (n_ok == 0L) return(1)  ## conservative: can't certify, treat as all-violating
      sum(mom_w$ok & (q_draws > mom_w$E_omega / mom_w$Var_omega)) / n_ok
    }

    w_lo <- 1e-6
    w_hi <- 0.5 - 1e-6
    pct_lo <- pct_outside_exact_for_w(w_lo)
    pct_hi <- pct_outside_exact_for_w(w_hi)

    if (pct_lo <= alpha_target) {
      ## Already below target with almost no pseudo-observations -- no
      ## sharpening needed at all.
      w_exact <- 0
      clipped_exact <- FALSE
    } else if (pct_hi > alpha_target) {
      ## Even the package's 0.5 per-group ceiling isn't enough.
      w_exact <- 0.5
      clipped_exact <- TRUE
    } else {
      g_exact <- function(w) pct_outside_exact_for_w(w) - alpha_target
      w_exact <- stats::uniroot(g_exact, lower = w_lo, upper = w_hi, tol = 1e-4)$root
      clipped_exact <- FALSE
    }
    w_exact_floored <- max(w_exact, w_floor)

    data.frame(
      group = lev, n_j = n_j,
      pct_outside_at_current = pct_outside_now,
      pct_outside_at_current_exact = pct_outside_now_exact,
      w_star_noncentral_untrunc = w_star, clipped_at_0.5_untrunc = clipped,
      w_star_floored_untrunc = w_star_floored,
      w_star_exact = w_exact, clipped_at_0.5_exact = clipped_exact,
      w_star_exact_floored = w_exact_floored,
      w_star_central_only = w_central
    )
  })
  do.call(rbind, rows)
}

.tmp_print_group_pwt_measurement_noncentral <- function(tab, alpha_target = 0.01, w_floor = 0.1) {
  cat(sprintf(
    "\n=== Group-specific pwt_measurement targeting alpha_target = %.3f (floor = %.2f) ===\n\n",
    alpha_target, w_floor
  ))
  print(tab[order(-tab$pct_outside_at_current_exact), ], row.names = FALSE, digits = 4)
  cat(sprintf(
    "\nGroups requiring clipping at the package's 0.5 per-group ceiling (exact criterion): %d / %d\n",
    sum(tab$clipped_at_0.5_exact), nrow(tab)
  ))
  cat(sprintf(
    "Groups where the floor (%.2f) actually binds (unfloored exact target below it): %d / %d\n",
    w_floor, sum(tab$w_star_exact < w_floor), nrow(tab)
  ))
  cat(sprintf(
    "Correlation(w_star_exact, w_star_noncentral_untrunc): %.3f\n",
    stats::cor(tab$w_star_exact, tab$w_star_noncentral_untrunc)
  ))
  cat(sprintf(
    "Correlation(w_star_exact, w_star_central_only): %.3f\n",
    stats::cor(tab$w_star_exact, tab$w_star_central_only)
  ))
  invisible(tab)
}

## ---------------------------------------------------------------------------
## Usage: with an existing ps/design already in your session (Ex_13/Ex_13b/
## Ex_13c's variable names: ps, design, grp, group_levels, re_names -- this
## is a PRE-RUN diagnostic, no 'fit' needed), just call:
##
##   tab_pwt <- .tmp_group_pwt_measurement_noncentral(
##     ps = ps, design = design, group = grp, group_levels = group_levels,
##     re_coef_names = re_names
##   )
##   .tmp_print_group_pwt_measurement_noncentral(tab_pwt)
##
##   ## Vector, ready to inspect (NOT fed into Prior_Setup_lmebayes() here):
##   w_vec <- stats::setNames(
##     tab_pwt$w_star_exact_floored[match(group_levels, tab_pwt$group)], group_levels
##   )
##   print(round(w_vec, 4))
##
## Columns:
##   pct_outside_at_current[_exact]     -- untruncated vs. exact pct outside
##                                          AT ps's CURRENT pwt_measurement.
##   w_star_noncentral_untrunc/         -- LEGACY: single-shot closed-form
##     clipped_at_0.5_untrunc/             inversion against the UNTRUNCATED
##     w_star_floored_untrunc              boundary (a0_j cancels out there).
##   w_star_exact/clipped_at_0.5_exact/ -- NEW: 1-D uniroot() over w against
##     w_star_exact_floored                the EXACT truncation-aware
##                                          boundary (Section 16.6).
##   w_star_central_only                -- no-shrinkage central-F baseline
##                                          (ignores noncentrality entirely).
##
## No demo re-run needed -- this file only *defines* the functions above;
## sourcing it does not execute anything else.
## ---------------------------------------------------------------------------
