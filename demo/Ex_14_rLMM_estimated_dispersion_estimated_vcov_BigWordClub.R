## Demo: rLMMindepNormalGamma_reg_estimated_vcov() called directly on
## bayesrules::big_word_club
##
## Case 5 of 5: ESTIMATED observation dispersion -- a *separate*, per-group
## sigma^2_j (dGamma() ING prior on each group's own precision, via
## dGamma_list()), NOT a single value shared across all groups -- AND
## ESTIMATED random-effect variance components (every Block~2 pfamily
## component is dIndependent_Normal_Gamma(), so tau^2_k is sampled too).
##
## Same model as demo("Ex_25_lmerb_dGamma_ING_BigWordClub", package =
## "lmebayes"), but this script calls rLMMindepNormalGamma_reg_estimated_vcov()
## directly instead of going through lmerb()/rlmerb(): model_setup(),
## Prior_Setup_GLMM(), pfamily_list(), and dGamma_list() (all exported
## from lmebayesCore) build the design and priors, then the script assembles
## by hand the exact 'group'/'prior_list' arguments that matrix_args_lmm()
## builds internally for rlmerb(), and calls the matrix-level export
## directly.
##
##   demo("Ex_14_rLMM_estimated_dispersion_estimated_vcov_BigWordClub", package = "lmebayesCore")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}

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

## ---------------------------------------------------------------------------
## 1. Per-group ING measurement dispersion requires every group to be full
##    column rank (no accept/reject envelope is built for rank-deficient
##    groups yet). Mirrors demo("Ex_25_lmerb_dGamma_ING_BigWordClub", package =
##    "lmebayes"): filter to full-rank schools, and drop school_id 2/18
##    (a known Block~1 per-group ING envelope sign-violation case, tracked
##    independently of this demo).
## ---------------------------------------------------------------------------
design_all <- model_setup(form_lmer, data = dat)
full_rank_schools <- names(design_all$groupef.rank)[design_all$groupef.rank]
cat(sprintf(
  "\n=== Full-rank filter: %d of %d schools kept ===\n",
  length(full_rank_schools),
  length(design_all$groupef.rank)
))
if (length(full_rank_schools) < length(design_all$groupef.rank)) {
  cat(
    "  Dropped:",
    paste(names(design_all$groupef.rank)[!design_all$groupef.rank], collapse = ", "),
    "\n"
  )
}
dat <- subset(dat, school_id %in% full_rank_schools)
dat$school_id <- droplevels(dat$school_id)

## TEMP: school 18 triggers ING envelope sign violation (UB2 < 0); drop and retest
temp_drop_schools <- c("18", "2")
drop <- intersect(temp_drop_schools, levels(dat$school_id))
if (length(drop)) {
  cat(sprintf(
    "\n=== TEMP: excluding school_id %s (Block~1 ING envelope failure) ===\n",
    paste(drop, collapse = ", ")
  ))
  dat <- subset(dat, !as.character(school_id) %in% drop)
  dat$school_id <- droplevels(dat$school_id)
}

## ---------------------------------------------------------------------------
## 2. Design + priors: model_setup() / Prior_Setup_GLMM() / pfamily_list()
##    / dGamma_list()
##
## pwt_dispersion = 0.2 calibrates the dIndependent_Normal_Gamma() Gamma
## window on each tau^2_k (wider/more diffuse than the pwt = 0.01 default).
## dispformula = ~school_id (matching the grouping factor exactly) requests
## Prior_Setup_GLMM()'s per-group Block~1 calibration
## (ing_prior_measurement_group), consumed below by dGamma_list().
## ---------------------------------------------------------------------------
design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup (full-rank schools only) ===\n\n")
print(design)
stopifnot(all(design$groupef.rank))

ps <- Prior_Setup_GLMM(
  form_lmer,
  data           = dat,
  pwt            = 0.01,
  pwt_dispersion = 0.2,
  dispformula    = ~school_id
)
cat("\n=== Prior_Setup_GLMM (ING + per-group Block~1 calibration) ===\n\n")
print(ps)

## Every Block~2 component is dIndependent_Normal_Gamma(): tau^2_k is
## *estimated* (sampled each sweep), not fixed at the lmer REML plug-in.
pf <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")

## One dGamma() pfamily per group level: each school_id gets its own
## sigma^2_j prior (shape/rate mean-matched to that school's own OLS/BLUP
## residual variance), not a shared pooled sigma^2.
disp_pf_list <- dGamma_list(ps)

## ---------------------------------------------------------------------------
## 3. Arguments matrix_args_lmm() would build for rlmerb() -- assembled here
##    by hand so rLMMindepNormalGamma_reg_estimated_vcov() can be called
##    directly.
##
## The routed export's 'prior_list' for a per-group ING Block~1 dispersion is
## NOT the dGamma() pfamily list itself -- it is a flat list with 'mu'/'Sigma'
## (the Block~2 hyperparameter prior, same shape .rLMM_validate_ing_
## measurement_prior_list() expects) plus 'shape_group'/'rate_group'/
## 'disp_lower_group'/'disp_upper_group' (one named-by-group-level numeric
## vector each), extracted here from each group's dGamma() pfamily --
## mirroring .lmebayes_resolve_group.dispersion_group_list() /
## .lmebayes_ing_measurement_prior_list_group() in mixed_rmerb_helpers.R.
## ---------------------------------------------------------------------------

## group_name is not a formal on the routed export; attach it to 'group'
## instead of relying on substitute() (see .lmebayes_resolve_group_name()).
grp <- design$group
attr(grp, "group_name") <- design$group_name

group_levels <- levels(grp)
re_names     <- design$groupef.names
p_re         <- length(re_names)

shape_group      <- stats::setNames(numeric(length(group_levels)), group_levels)
rate_group       <- stats::setNames(numeric(length(group_levels)), group_levels)
disp_lower_group <- stats::setNames(numeric(length(group_levels)), group_levels)
disp_upper_group <- stats::setNames(numeric(length(group_levels)), group_levels)
for (lev in group_levels) {
  pl <- disp_pf_list[[lev]]$prior_list
  shape_group[[lev]]      <- pl$shape[1L]
  rate_group[[lev]]       <- pl$rate[1L]
  disp_lower_group[[lev]] <- pl$disp_lower
  disp_upper_group[[lev]] <- pl$disp_upper
}

prior_list <- list(
  mu               = matrix(0, nrow = p_re, ncol = 1L, dimnames = list(re_names, NULL)),
  Sigma            = as.matrix(ps$group.Sigma),
  shape_group      = shape_group,
  rate_group       = rate_group,
  disp_lower_group = disp_lower_group,
  disp_upper_group = disp_upper_group
)

cat(sprintf(
  "\n=== Per-group sigma^2_j ING prior: %d groups, sigma^2_hat range [%.4f, %.4f] ===\n\n",
  length(group_levels),
  min(rate_group / (shape_group - 1)),
  max(rate_group / (shape_group - 1))
))

## ---------------------------------------------------------------------------
## 4. lmer reference fit
## ---------------------------------------------------------------------------
cat("\n=== lmer reference fit ===\n\n")
fit_lmer <- lme4::lmer(form_lmer, data = dat, REML = TRUE)
print(summary(fit_lmer))
cat(sprintf(
  "\n  Pooled REML sigma^2 (compare to per-group sigma^2_j below): %.4f\n",
  stats::sigma(fit_lmer)^2
))

## ---------------------------------------------------------------------------
## 5. Direct call: rLMMindepNormalGamma_reg_estimated_vcov()
##
## Two-block Gibbs with both an ING Block~1 sweep (sigma^2_j estimated) and
## an ING Block~2 sweep (tau^2_k estimated): an optional pilot stage
## (gap_tol/mode_gap_max) recenters the main stage's starting point away
## from the ICM mode, then Theorem~3 calibrates the number of inner sweeps
## per stored draw.
##
## progbar/verbose match demo("Ex_25_lmerb_dGamma_ING_BigWordClub", package =
## "lmebayes"): that demo calls lmerb() without overriding progbar/verbose,
## and lmerb()'s own formals are progbar = NULL (falsy -- no bar shown) and
## a hardcoded verbose = TRUE passed to rlmerb().
## ---------------------------------------------------------------------------
fit <- rLMMindepNormalGamma_reg_estimated_vcov(
  n            = 3000L,
  y            = design$y,
  D            = design$D,
  group        = grp,
  W            = design$W,
  pfamily_list = pf,
  dispprior_list = prior_list,
  gap_tol      = 0.05,
  mode_gap_max = 1.0,
  diag_sweeps  = FALSE,
  progbar      = FALSE,
  verbose      = TRUE
)

stopifnot(isTRUE(fit$any_non_normal))
stopifnot(is.matrix(fit$group.dispersion))
stopifnot(all(is.finite(fit$group.dispersion)), all(fit$group.dispersion > 0))
stopifnot(!is.null(fit$pilot$chisq))
stopifnot(!is.null(fit$pilot) && !is.null(fit$pilot$draws$sweep_history))
stopifnot(!is.null(fit$sweep_history))

n_draws <- nrow(fit$popef[[re_names[1L]]])

cat(sprintf(
  "\nPilot vs ICM mode (chi-squared): p = %.4g (n_pilot = %d, m_convergence_pilot = %s)\n",
  fit$pilot$chisq$p_value,
  fit$pilot$chisq$n_pilot,
  format(fit$pilot$m_convergence)
))
cat(sprintf("m_convergence (main) = %d\n", fit$m_convergence))

## ---------------------------------------------------------------------------
## 6. sigma^2_j (Block~1): post-sweep draws' per-group means vs
##    OLS/BLUP-calibrated prior mean, and vs the pooled lmer REML sigma^2.
## ---------------------------------------------------------------------------
## Build every row first (no printing) so the header/rows below print as one
## contiguous block -- otherwise demo()'s echo = TRUE interleaves this loop's
## own source between the header and the first data row.
disp_prior_mean <- rate_group / (shape_group - 1)
rows_disp <- character(0L)
for (lev in group_levels) {
  rows_disp <- c(rows_disp, sprintf(
    "  %-6s  %10.4f  %10.4f  [%.4f, %.4f]\n",
    lev, fit$group.dispersion.mean[[lev]], disp_prior_mean[[lev]],
    disp_lower_group[[lev]], disp_upper_group[[lev]]
  ))
}

## Single cat() call for the whole block (title + header + rows): demo()'s
## echo = TRUE prints one "> ..." per top-level R statement, so splitting the
## table across multiple cat() calls -- even adjacent ones -- still
## interleaves an echoed statement between each piece. One call in, one
## block out.
cat(
  "\n=== sigma^2_j: post mean vs calibrated prior mean (all groups) ===\n\n",
  sprintf("  %-6s  %10s  %10s  %10s\n",
          "group", "post mean", "prior mean", "[window]"),
  rows_disp,
  sep = ""
)

## ---------------------------------------------------------------------------
## 7. tau^2_k (Block~2): post-sweep draws stayed inside the calibrated Gamma
##    window.
## ---------------------------------------------------------------------------
cat("\n=== tau^2_k: post mean vs calibrated window ===\n\n")
for (k in re_names) {
  pr_k <- pf[[k]]$prior_list
  t2   <- fit$popef.dispersion[, k]
  cat(sprintf(
    "  %-18s post mean = %8.4f  [window (%.4f, %.4f)]\n",
    k, mean(t2), pr_k$disp_lower, pr_k$disp_upper
  ))
}

## ---------------------------------------------------------------------------
## 8. Block 2 fixed effects: Gibbs (MCMC) vs ICM mean vs glmmTMB fixef
##    (+ uncertainty)
##
## No exact-iid engine exists for this model (both sigma^2_j and tau^2_k are
## estimated, so the joint posterior is not exactly Gaussian) -- 'gibbs
## mean'/'gibbs SD' (the lmebayesCore output: posterior mean/SD of gamma_k
## across the main-stage MCMC draws) is compared directly to the same
## dispformula = ~school_id glmmTMB reference fit that calibrated
## Prior_Setup_GLMM()'s per-group Block~1 prior (ps$fit_ref), same
## column layout as Ex_11's Section 6 (just without the 'iid' columns, since
## no iid engine exists here). 'diff(SE)' re-expresses the gibbs-mean vs
## glmmTMB gap in units of glmmTMB's own Std. Error -- the right scale to
## judge it on.
## ---------------------------------------------------------------------------
fit_ref <- ps$fit_ref
fe_ref  <- lmebayesCore:::.lmebayes_reference_fixef(fit_ref)
se_ref  <- sqrt(diag(lmebayesCore:::.lmebayes_reference_vcov(fit_ref)))

## Build every row first (no printing) so the header/dashes/rows below print
## as one contiguous block -- otherwise demo()'s echo = TRUE interleaves this
## loop's own source between the header and the first data row.
rows_fe <- character(0L)
for (k in re_names) {
  dm_gibbs <- colMeans(fit$popef[[k]])
  sd_gibbs <- apply(fit$popef[[k]], 2L, sd)
  icm_k    <- fit$popef.mode[[k]]
  for (nm in names(dm_gibbs)) {
    fe_nm <- if (identical(k, "(Intercept)") && identical(nm, "(Intercept)")) {
      "(Intercept)"
    } else if (identical(nm, "(Intercept)")) {
      k
    } else if (identical(k, "(Intercept)")) {
      nm
    } else {
      cand <- c(paste0(nm, ":", k), paste0(k, ":", nm))
      hit  <- cand[cand %in% names(fe_ref)]
      if (length(hit)) hit[1L] else NA_character_
    }
    fe_val <- if (!is.na(fe_nm) && fe_nm %in% names(fe_ref)) unname(fe_ref[fe_nm]) else NA_real_
    se_val <- if (!is.na(fe_nm) && fe_nm %in% names(se_ref)) unname(se_ref[fe_nm]) else NA_real_
    diff_se <- (dm_gibbs[[nm]] - fe_val) / se_val
    rows_fe <- c(rows_fe, sprintf(
      "  %-18s  %-28s  %10.4f  %8.4f  %10.4f  %10.4f  %8.4f  %8.2f\n",
      k, nm, dm_gibbs[[nm]], sd_gibbs[[nm]], icm_k[[nm]], fe_val, se_val, diff_se
    ))
  }
}

## Single cat() call for the whole block (title + header + dashes + rows +
## footnote): demo()'s echo = TRUE prints one "> ..." per top-level R
## statement, so splitting the table across multiple cat() calls -- even
## adjacent ones -- still interleaves an echoed statement between each
## piece. One call in, one block out.
cat(
  "\n=== Block 2 fixed effects: Gibbs vs ICM mean vs glmmTMB fixef (+ uncertainty) ===\n\n",
  sprintf("  %-18s  %-28s  %10s  %8s  %10s  %10s  %8s  %8s\n",
          "RE component", "parameter", "gibbs mean", "gibbs SD",
          "ICM mean", "glmmTMB", "glmm SE", "diff(SE)"),
  sprintf("  %-18s  %-28s  %10s  %8s  %10s  %10s  %8s  %8s\n",
          strrep("-", 18L), strrep("-", 28L), strrep("-", 10L), strrep("-", 8L),
          strrep("-", 10L), strrep("-", 10L), strrep("-", 8L), strrep("-", 8L)),
  rows_fe,
  "\n  diff(SE) = (gibbs mean - glmmTMB estimate) / glmmTMB Std. Error -- |diff(SE)| < ~1-2\n",
  "  is well within glmmTMB's own uncertainty for that coefficient, not a discrepancy.\n",
  sep = ""
)

## ---------------------------------------------------------------------------
## 8b. Extended (lambda + Omega)-aware local rate diagnostic
##     (inst/BLOCK_GIBBS_ERGODICITY_ING.md Sections 6, 14, 15)
##
## Both extensions of the theory note apply here (tau^2_k *and* sigma^2_j
## are estimated): lambda_ing (Section 6) and omega_ing (Section 14) are
## supplied together; Sections 14.5/15 show the two stack with no new
## cross-derivation (P11^(gamma,Omega) = P11^(lambda,Omega) = 0 exactly).
##
## lambda_ing/omega_ing's reference precisions and residuals are ALL
## evaluated at the now-*completed* sampler's own posterior-mean state
## (fit$popef.dispersion for tau^2_k, fit$group.dispersion.mean for
## sigma^2_j, fit$groupef/fit$popef for beta_j/gamma_k), rather than
## the certified corner's P_b/disp_upper_group plug-ins or an external
## dispformula = ~school_id glmmTMB reference fit (fit_ref). Every quantity
## below therefore comes from the SAME fitted object and describes one
## mutually-consistent reference state; mixing corner/external-reference-fit
## quantities can put a precision and its paired residual badly out of step
## for a group/component the reference poorly predicts -- the joint
## (beta, lambda/Omega) Hessian is only guaranteed PD near its own mode, see
## two_block_rate_ing()'s "local, uncertified diagnostic" documentation and
## demo("Ex_13_...")'s Section 7b note for the worked-out det(H) < 0 case
## this fix followed from. This is a *diagnostic*, not a certified bound:
## see two_block_rate_ing()'s own documentation for why no analogue of
## Theorem~3/Corollary~1 exists for the extended chain.
## ---------------------------------------------------------------------------
lambda_post <- stats::setNames(
  1 / vapply(re_names, function(k) mean(fit$popef.dispersion[, k]), numeric(1L)),
  re_names
)
u_post <- lmebayesCore:::.lmebayes_posterior_u(
  fit, design$group_name, re_names, design$W
)

lambda_ing <- stats::setNames(lapply(re_names, function(k) {
  list(
    lambda = unname(lambda_post[[k]]),
    shape  = pf[[k]]$prior_list$shape[[1L]],
    u      = stats::setNames(u_post[[k]], rownames(u_post))
  )
}), re_names)

n_group  <- stats::setNames(as.numeric(table(grp)), group_levels)
e_post   <- lmebayesCore:::.lmebayes_posterior_group_residuals(
  fit, y = design$y, D = design$D, group = grp,
  group_name = design$group_name, groupef.names = re_names
)
omega_post <- stats::setNames(
  1 / fit$group.dispersion.mean[group_levels], group_levels
)
omega_ing <- list(
  scope = "per_group",
  omega = omega_post,
  shape = shape_group[group_levels],
  n     = n_group,
  e     = e_post
)

## P_b_ref (the calibrated Block~1 prior precision the certified corner
## rate/m_convergence calibration itself uses) is kept as-is for
## prior_list_block1_rate$P -- unlike lambda_ing/omega_ing's own reference
## precisions above, this is NOT switched to the posterior mean, so
## lambda_star_base below still reproduces that same certified corner rate
## (mirrors demo("Ex_13_...")'s Section 7b keeping disp_upper_group for the
## same reason).
P_b_ref <- lmebayesCore:::.rLMM_P_from_pfamily_list(pf, re_names)

prior_list_block1_rate <- list(P = P_b_ref, dispersion = disp_upper_group[group_levels])
prior_list_block2_rate <- lapply(pf, function(pfk) {
  pl <- pfk$prior_list
  list(
    mu = pl$mu, Sigma = pl$Sigma,
    dispersion = if (identical(pfk$pfamily, "dNormal")) pl$dispersion else pl$disp_lower
  )
})

rate_ext <- two_block_rate_ing(
  x = design$D, group = grp, x_hyper = design$W,
  prior_list_block1 = prior_list_block1_rate,
  prior_list_block2 = prior_list_block2_rate,
  lambda_ing = lambda_ing,
  omega_ing  = omega_ing
)
cat("\n=== Extended (lambda + Omega)-aware local rate diagnostic (Sections 6, 14-15) ===\n\n")
print(rate_ext)

## ---------------------------------------------------------------------------
## 8c. Empirical worst-case: same (lambda + Omega)-aware extended rate,
##     evaluated at EVERY main-stage draw instead of the single
##     posterior-mean reference state above.
##
## Both tau^2_k and sigma^2_j are estimated here, so each of the n
## main-stage draws has its own sampled tau^2_k (fit$popef.dispersion[i, ])
## AND sigma^2_j (fit$group.dispersion[i, ]), together with its own
## beta_j/gamma_p (fit$groupef/fit$popef draw i) -- and hence its own
## u_jp/e_j. Treating the n draws as approximate posterior samples, this
## asks: what is the largest local rate actually realized across them,
## rather than at one plug-in point? Mirrors the pilot-draw pmax() scan
## (.two_block_pilot_ub_from_coefficients()) but for the *extended* system
## and the *main*-stage output; still a diagnostic, not a certified bound
## (inst/BLOCK_GIBBS_ERGODICITY_ING.md Sections 7, 9, 11, 15).
## ---------------------------------------------------------------------------
lambda_spec <- stats::setNames(lapply(re_names, function(k) {
  list(shape = pf[[k]]$prior_list$shape[[1L]])
}), re_names)
omega_spec <- list(
  scope = "per_group",
  shape = shape_group[group_levels],
  n     = n_group
)

rate_emp <- lmebayesCore:::.two_block_rate_ing_over_draws(
  fit = fit, n_draws = n_draws,
  x = design$D, group = grp, x_hyper = design$W,
  prior_list_block1 = prior_list_block1_rate,
  prior_list_block2 = prior_list_block2_rate,
  group_name = design$group_name, groupef.names = re_names,
  y = design$y, D = design$D,
  lambda_spec = lambda_spec, omega_spec = omega_spec
)

cat(sprintf(
  paste0(
    "\n=== Empirical (lambda + Omega)-aware local rate over n = %d main-stage draws ===\n\n",
    "  lambda* (extended, single reference state)     = %.6f\n",
    "  lambda* (extended, empirical max over %d draws) = %.6f  (draw #%d)\n",
    "  draws with lambda*(extended) >= 1: %d / %d\n"
  ),
  rate_emp$n_draws,
  rate_ext$lambda_star,
  rate_emp$n_draws, rate_emp$lambda_star_max, rate_emp$i_max,
  rate_emp$n_over_one, rate_emp$n_draws
))

## Combined mean-bias/Var_final-ratio charts (Claims 1 and 3 of the two-block
## Gibbs ergodicity reference): rLMMindepNormalGamma_reg_estimated_vcov() goes
## through the sweeps-outer/chains-inner pilot/main engine, so both
## fit$pilot$draws$sweep_history and fit$sweep_history carry cov_by_sweep (as does
## rLMMNormal_reg_known_vcov(sim_method = "TWO_BLOCK_GIBBS")'s single main
## stage now that its engine is rGLMM_sweep()-based too). Both per-group
## dispersion and group.Sigma are *estimated* here (that is the point of this
## demo), so plot_mean_convergence()/plot_var_convergence()'s fit-object
## methods automatically find no exact reference available and fall back to
## the empirical mean_final/Var_final (last-sweep cross-chain mean/
## covariance). 'stage' selects fit$pilot$draws$sweep_history ("pilot") or
## fit$sweep_history ("main"); n_chains defaults to the correct chain count
## for each stage (fit$pilot$chisq$n_pilot vs fit$n). The pilot stage's whole
## job is to *recenter* the main stage's starting point away from the ICM
## mode, so the mean-bias chart is the more direct diagnostic there; the main
## stage is where the variance/uncertainty calibration check matters most --
## both charts are shown for both stages for consistency.
for (stg in c("pilot", "main")) {
  plot_mean_convergence(fit, whitened = FALSE, stage = stg)
  plot_mean_convergence(fit, whitened = TRUE, stage = stg)
  plot_var_convergence(fit, whitened = FALSE, stage = stg)
  plot_var_convergence(fit, whitened = TRUE, stage = stg)
}

## component = "precision" plots the RE variance-covariance itself,
## group.Sigma = diag(tau2_k), tracked as its reciprocal 1/tau2_k
## (precision) rather than tau2_k -- E[1/tau2_k] stays well-defined under a
## weak/vague prior where E[tau2_k] need not be. Unlike the Block~2
## fixed-effects charts above, there is no exact reference here
## (design/measurement_prior_list are ignored for this component) and no
## whitened mode yet (only each component's own cross-chain variance is
## captured so far, not its covariance with other components/fixed effects
## -- see inst/BLOCK_GIBBS_ERGODICITY.md's "Future work" section), so only
## the named (non-whitened) charts are drawn.
for (stg in c("pilot", "main")) {
  plot_mean_convergence(fit, component = "precision", stage = stg)
  plot_var_convergence(fit, component = "precision", stage = stg)
}

## ---------------------------------------------------------------------------
## 9. Block 2 hyperparameters: prior mean, ICM (gamma @ lmer tau2), pilot
##    mean, MCMC mean (supplementary to Section 8 -- shows how far the
##    pilot-stage recentering moved the starting point away from the prior).
## ---------------------------------------------------------------------------
cn <- unlist(lapply(re_names, function(k) {
  paste0(k, "::", colnames(fit$popef[[k]]))
}))
beta_bar    <- unlist(lapply(re_names, function(k) colMeans(fit$popef[[k]])))
theta_icm   <- unlist(lapply(re_names, function(k) fit$popef.mode[[k]]))
theta_prior <- unlist(lapply(re_names, function(k) {
  nms <- colnames(fit$popef[[k]])
  ## Raw pfamily_list() objects (unlike lmerb()'s processed fit$prior) store
  ## the Block~2 prior mean as prior_list$mu, an ncol(W[[k]]) x 1 matrix
  ## dimnamed by colnames(W[[k]]) -- not prior_list$mu.
  unname(pf[[k]]$prior_list$mu[nms, 1L])
}))
theta_pilot <- unlist(lapply(re_names, function(k) {
  nms <- colnames(fit$popef[[k]])
  unname(fit$popef.init[[k]][nms])
}))
names(beta_bar) <- names(theta_icm) <- names(theta_prior) <- names(theta_pilot) <- cn

block2_cmp <- data.frame(
  prior_mean      = unname(theta_prior),
  icm_lmer_tau2   = unname(theta_icm),
  pilot_mean      = unname(theta_pilot),
  mcmc_mean       = unname(beta_bar),
  row.names       = cn,
  check.names     = FALSE
)
cat("\n=== Block 2 hyperparameters (prior / ICM / pilot / MCMC) ===\n\n")
print(round(block2_cmp, 4))

## ---------------------------------------------------------------------------
## 10. Sweep-history diagnostics: cross-chain mean/SD vs inner sweep, for
##    both the pilot and main stages.
##
## plot_sweep_history_diag(engine = "base") stacks one panel per coef_focus
## entry via mfrow = c(length(coef_focus), 1L) -- passing all 7 at once needs
## a very tall plot device and can error with "figure margins too large" on
## an ordinary-sized device. Mirrors demo("Ex_19_glmerb_book_banning_state_
## covariates", package = "lmebayes"): split coef_focus into small batches
## (2-3 panels) and call plot_sweep_history_diag() once per batch.
## ---------------------------------------------------------------------------
coef_focus_all <- list(
  c("(Intercept)", "(Intercept)"),
  c("(Intercept)", "private_school"),
  c("(Intercept)", "title1"),
  c("(Intercept)", "free_reduced_lunch"),
  c("distracted_ppvt", "(Intercept)"),
  c("distracted_a1", "(Intercept)"),
  c("distracted_a1", "free_reduced_lunch")
)
coef_focus_batches <- list(
  coef_focus_all[1:2],
  coef_focus_all[3:4],
  coef_focus_all[5:7]
)

for (stg in c("pilot", "main")) {
  for (batch in coef_focus_batches) {
    plot_sweep_history_diag(fit, batch, stage = stg)
  }
}

## ---------------------------------------------------------------------------
## 11. Random effects: MCMC mean (per group, per draw average) vs ICM mode
##
## fit$groupef: long data.frame with one row per (draw, group) -- beta_j
## draws (the full, non-centered coefficient; see ?rLMM_reg's "Model and
## notation" section). fit$groupef.mode is the ICM mode these Gibbs sweeps
## started from.
## ---------------------------------------------------------------------------
grp_col  <- design$group_name
grp_levs <- rownames(fit$groupef.mode)

re_draws_mean <- tapply(
  seq_len(nrow(fit$groupef)),
  fit$groupef[[grp_col]],
  function(idx) colMeans(fit$groupef[idx, re_names, drop = FALSE]),
  simplify = FALSE
)
re_draws_sd <- tapply(
  seq_len(nrow(fit$groupef)),
  fit$groupef[[grp_col]],
  function(idx) apply(fit$groupef[idx, re_names, drop = FALSE], 2L, sd),
  simplify = FALSE
)

## Build every row first (no printing) so the header/dashes/rows below print
## as one contiguous block -- otherwise demo()'s echo = TRUE interleaves this
## loop's own source between the header and the first data row.
rows_re <- character(0L)
for (lev in grp_levs) {
  lev_chr <- as.character(lev)
  for (k in re_names) {
    mcmc_m <- re_draws_mean[[lev_chr]][[k]]
    mcmc_s <- re_draws_sd[[lev_chr]][[k]]
    icm_m  <- fit$groupef.mode[lev_chr, k]
    se_val <- mcmc_s / sqrt(n_draws)
    z_val  <- (mcmc_m - icm_m) / se_val
    rows_re <- c(rows_re, sprintf(
      "  %-6s  %-18s  %10.4f  %10.4f  %10.4f  %6.2f\n",
      lev_chr, k, mcmc_m, icm_m, se_val, z_val
    ))
  }
}

## Single cat() call for the whole block (title + header + dashes + rows +
## footnote): demo()'s echo = TRUE prints one "> ..." per top-level R
## statement, so splitting the table across multiple cat() calls -- even
## adjacent ones -- still interleaves an echoed statement between each
## piece. One call in, one block out.
cat(
  "\n=== Random effects: MCMC mean vs ICM mode (all groups) ===\n\n",
  sprintf("  %-6s  %-18s  %10s  %10s  %10s  %6s\n",
          "group", "RE component", "MCMC mean", "ICM mode", "SE(mean)", "z"),
  sprintf("  %-6s  %-18s  %10s  %10s  %10s  %6s\n",
          strrep("-", 6L), strrep("-", 18L),
          strrep("-", 10L), strrep("-", 10L), strrep("-", 10L), strrep("-", 6L)),
  rows_re,
  "\n  Note: MCMC draws here are autocorrelated (two-block Gibbs, not the\n",
  "  exact-iid engine), so z far from 0 does not by itself indicate a\n",
  "  problem -- treat these as approximate; see convergence_info for the\n",
  "  Theorem~3 inner-sweep count used.\n",
  sep = ""
)
