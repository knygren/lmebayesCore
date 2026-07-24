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
## Prior_Setup_lmebayes(), pfamily_list(), and dGamma_list() (all exported
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
full_rank_schools <- names(design_all$re_rank)[design_all$re_rank]
cat(sprintf(
  "\n=== Full-rank filter: %d of %d schools kept ===\n",
  length(full_rank_schools),
  length(design_all$re_rank)
))
if (length(full_rank_schools) < length(design_all$re_rank)) {
  cat(
    "  Dropped:",
    paste(names(design_all$re_rank)[!design_all$re_rank], collapse = ", "),
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
## 2. Design + priors: model_setup() / Prior_Setup_lmebayes() / pfamily_list()
##    / dGamma_list()
##
## pwt_dispersion = 0.2 calibrates the dIndependent_Normal_Gamma() Gamma
## window on each tau^2_k (wider/more diffuse than the pwt = 0.01 default).
## dispformula = ~school_id (matching the grouping factor exactly) requests
## Prior_Setup_lmebayes()'s per-group Block~1 calibration
## (ing_prior_measurement_group), consumed below by dGamma_list().
## ---------------------------------------------------------------------------
design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup (full-rank schools only) ===\n\n")
print(design)
stopifnot(all(design$re_rank))

ps <- Prior_Setup_lmebayes(
  form_lmer,
  data           = dat,
  pwt            = 0.01,
  pwt_dispersion = 0.2,
  dispformula    = ~school_id
)
cat("\n=== Prior_Setup_lmebayes (ING + per-group Block~1 calibration) ===\n\n")
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
## mirroring .lmebayes_resolve_dispersion_ranef_group_list() /
## .lmebayes_ing_measurement_prior_list_group() in mixed_rmerb_helpers.R.
## ---------------------------------------------------------------------------

## group_name is not a formal on the routed export; attach it to 'group'
## instead of relying on substitute() (see .lmebayes_resolve_group_name()).
grp <- design$groups
attr(grp, "group_name") <- design$group_name

group_levels <- levels(grp)
re_names     <- design$re_coef_names
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
  Sigma            = as.matrix(ps$Sigma_ranef),
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
  D            = design$Z,
  group        = grp,
  W            = design$X_hyper,
  prior_list   = prior_list,
  pfamily_list = pf,
  gap_tol      = 0.05,
  mode_gap_max = 1.0,
  diag_sweeps  = FALSE,
  progbar      = FALSE,
  verbose      = TRUE
)

stopifnot(isTRUE(fit$any_non_normal))
stopifnot(is.matrix(fit$dispersion_ranef))
stopifnot(all(is.finite(fit$dispersion_ranef)), all(fit$dispersion_ranef > 0))
stopifnot(!is.null(fit$pilot_chisq))
stopifnot(!is.null(fit$pilot) && !is.null(fit$pilot$sweep_history))
stopifnot(!is.null(fit$sweep_history))

n_draws <- nrow(fit$fixef[[re_names[1L]]])

cat(sprintf(
  "\nPilot vs ICM mode (chi-squared): p = %.4g (n_pilot = %d, m_convergence_pilot = %s)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot,
  format(fit$m_convergence_pilot)
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
    lev, fit$dispersion_ranef.mean[[lev]], disp_prior_mean[[lev]],
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
  t2   <- fit$fixef.dispersion[, k]
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
## Prior_Setup_lmebayes()'s per-group Block~1 prior (ps$fit_ref), same
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
  dm_gibbs <- colMeans(fit$fixef[[k]])
  sd_gibbs <- apply(fit$fixef[[k]], 2L, sd)
  icm_k    <- fit$fixef.mode[[k]]
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

## Combined Var/Var_final ratio chart (Claim 3 of the two-block Gibbs
## ergodicity reference): rLMMindepNormalGamma_reg_estimated_vcov() goes
## through the sweeps-outer/chains-inner pilot/main engine, so both
## fit$pilot$sweep_history and fit$sweep_history carry cov_by_sweep -- unlike
## rLMMNormal_reg_known_vcov(sim_method = "TWO_BLOCK_GIBBS"), which does not
## capture sweep history yet. Both per-group dispersion and Sigma_ranef are
## *estimated* here (that is the point of this demo), so no exact reference
## covariance is available: 'design'/'measurement_prior_list' are
## intentionally omitted below and the plot falls back to the empirical
## Var_final (last-sweep cross-chain covariance).
for (st in list(fit$pilot$sweep_history, fit$sweep_history)) {
  if (is.null(st)) next
  plot_sweep_history_var_ratio(st, whitened = FALSE)
  plot_sweep_history_var_ratio(st, whitened = TRUE)
}

## ---------------------------------------------------------------------------
## 9. Block 2 hyperparameters: prior mean, ICM (gamma @ lmer tau2), pilot
##    mean, MCMC mean (supplementary to Section 8 -- shows how far the
##    pilot-stage recentering moved the starting point away from the prior).
## ---------------------------------------------------------------------------
cn <- unlist(lapply(re_names, function(k) {
  paste0(k, "::", colnames(fit$fixef[[k]]))
}))
beta_bar    <- unlist(lapply(re_names, function(k) colMeans(fit$fixef[[k]])))
theta_icm   <- unlist(lapply(re_names, function(k) fit$fixef.mode[[k]]))
theta_prior <- unlist(lapply(re_names, function(k) {
  nms <- colnames(fit$fixef[[k]])
  ## Raw pfamily_list() objects (unlike lmerb()'s processed fit$prior) store
  ## the Block~2 prior mean as prior_list$mu, an ncol(W[[k]]) x 1 matrix
  ## dimnamed by colnames(W[[k]]) -- not prior_list$mu_fixef.
  unname(pf[[k]]$prior_list$mu[nms, 1L])
}))
theta_pilot <- unlist(lapply(re_names, function(k) {
  nms <- colnames(fit$fixef[[k]])
  unname(fit$fixef.init[[k]][nms])
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

for (st in list(fit$pilot$sweep_history, fit$sweep_history)) {
  if (is.null(st)) next
  for (batch in coef_focus_batches) {
    plot_sweep_history_diag(st, batch)
  }
}

## ---------------------------------------------------------------------------
## 11. Random effects: MCMC mean (per group, per draw average) vs ICM mode
##
## fit$coefficients: long data.frame with one row per (draw, group) -- beta_j
## draws (the full, non-centered coefficient; see ?rLMM_reg's "Model and
## notation" section). fit$ranef.mode is the ICM mode these Gibbs sweeps
## started from.
## ---------------------------------------------------------------------------
grp_col  <- design$group_name
grp_levs <- rownames(fit$ranef.mode)

re_draws_mean <- tapply(
  seq_len(nrow(fit$coefficients)),
  fit$coefficients[[grp_col]],
  function(idx) colMeans(fit$coefficients[idx, re_names, drop = FALSE]),
  simplify = FALSE
)
re_draws_sd <- tapply(
  seq_len(nrow(fit$coefficients)),
  fit$coefficients[[grp_col]],
  function(idx) apply(fit$coefficients[idx, re_names, drop = FALSE], 2L, sd),
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
    icm_m  <- fit$ranef.mode[lev_chr, k]
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
