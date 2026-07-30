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
## Key simplification (Section 9.3): sigma2_hat_j, hence Omega_hat_j, B_j,
## beta_bar_j, and the whole DISTRIBUTION of q_j = (beta_j-beta_hat_ols_j)'
## D_j'D_j(beta_j-beta_hat_ols_j) is invariant to pwt_measurement -- only the
## boundary K_j = 2*r0_j+RSS_ols_j moves as pwt_measurement (hence r0_j)
## changes. So ONE Monte Carlo simulation of q_j's distribution per group is
## enough to solve for w_j via a single quantile lookup + closed-form
## algebra, with no need to re-simulate per candidate w_j and no need for
## iterative root-finding.
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

    ing_j       <- ps$ing_prior_measurement_group[[lev]]
    sigma2_hat  <- ing_j$sigma2_hat
    Omega_hat_j <- 1 / sigma2_hat
    rate_now    <- ing_j$rate           ## r0_j at ps's CURRENT pwt_measurement
    thresh_now  <- 2 * rate_now + RSS_ols

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

    pct_outside_now <- 100 * mean(q_draws > thresh_now)

    ## Quantile-match to alpha_target, then invert r0 -> a0 -> n_prior -> w.
    q_star   <- as.numeric(stats::quantile(q_draws, probs = 1 - alpha_target, type = 7))
    r0_star  <- (q_star - RSS_ols) / 2
    a0_star  <- r0_star / sigma2_hat + 1
    n_prior_star <- 2 * a0_star - p_re - 1
    w_raw    <- n_prior_star / (n_prior_star + n_j)
    w_star   <- max(0, min(0.5, w_raw))
    clipped  <- w_raw > 0.5
    w_star_floored <- max(w_star, w_floor)

    w_central <- max(0, (nu_star_central - n_j - 1) / (nu_star_central - 1))

    data.frame(
      group = lev, n_j = n_j,
      pct_outside_at_current = pct_outside_now,
      w_star_noncentral = w_star, clipped_at_0.5 = clipped,
      w_star_floored = w_star_floored,
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
  print(tab[order(-tab$pct_outside_at_current), ], row.names = FALSE, digits = 4)
  cat(sprintf(
    "\nGroups requiring clipping at the package's 0.5 per-group ceiling: %d / %d\n",
    sum(tab$clipped_at_0.5), nrow(tab)
  ))
  cat(sprintf(
    "Groups where the floor (%.2f) actually binds (unfloored target below it): %d / %d\n",
    w_floor, sum(tab$w_star_noncentral < w_floor), nrow(tab)
  ))
  cat(sprintf(
    "Correlation(w_star_noncentral, w_star_central_only): %.3f\n",
    stats::cor(tab$w_star_noncentral, tab$w_star_central_only)
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
##     tab_pwt$w_star_floored[match(group_levels, tab_pwt$group)], group_levels
##   )
##   print(round(w_vec, 4))
##
## No demo re-run needed -- this file only *defines* the functions above;
## sourcing it does not execute anything else.
## ---------------------------------------------------------------------------
