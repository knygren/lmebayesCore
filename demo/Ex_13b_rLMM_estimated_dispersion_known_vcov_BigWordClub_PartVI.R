## Demo: rLMMindepNormalGamma_reg_known_vcov() called directly on
## bayesrules::big_word_club -- Part VI (model-derived) alternate prior
##
## Variant of demo("Ex_13_rLMM_estimated_dispersion_known_vcov_BigWordClub",
## package = "lmebayesCore"): SAME model, data, filtering, and sampler call --
## only the per-group Block~1 measurement-dispersion prior differs. Ex_13
## uses dGamma_list(ps, disp_upper_anchor = "symmetric") directly (Part
## I/II of inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md: mu_j held fixed at
## group j's own null-model intercept). This demo instead applies the
## un-implemented Part VI extension of that document -- "also integrating
## out the prior mean mu_j" -- with a *model-derived* Omega_j (not an
## invented one): mu_j[k] stands in for the model's own conditional prior
## mean for b_j, W_j[[k]] %*% gamma_k (the Block~2 hyper-regression).
## gamma_k's own prior uncertainty is already calibrated at
## ps$prior_list[[k]]$Sigma_fixef (the same object pfamily_list(ps) uses to
## build pf below), so the "fixed, known uncertainty about mu_j" Part VI
## asks for is exactly that existing Sigma_fixef propagated through group
## j's own hyper-design row:
##
##   Omega_j[k,k] = W_j[[k]] %*% Sigma_fixef_k %*% t(W_j[[k]])   (diagonal
##                                                      across components k)
##
## folded in as Sigma_j' = Sigma_j + Omega_j before the same
## compute_gaussian_prior() calibration Part I/dGamma_list() already use,
## then the same symmetric-quantile window construction
## (shape_w/rate_w from n_combined,j) applied to the resulting
## sigma2_hat,j'.
##
## UPDATE: Part VI's Omega_j fold-in is now wired into
## Prior_Setup_lmebayes() itself (permanent default, not opt-in), and
## disp_lower/disp_upper are now literal quantiles of the resulting
## Gamma(shape_ING, rate) (not the shape_w/rate_w/n_combined window below) --
## see inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md Part VI and R/dGamma_list_
## lmebayes_prior_setup.R. A plain `dGamma_list(ps)` call now reproduces
## the Omega_j fold-in automatically; this demo's hand-rolled block below
## (part_vi_group, via the internal helpers
## .lmebayes_ing_prior_measurement_group_glm_inputs() and
## .lmebayes_compute_ing_prior_cal_from_sigma()) is kept as a from-scratch
## worked derivation for comparison, and its shape/rate should now match
## ps$ing_prior_measurement_group's own values (only the window construction
## still differs -- see the comparison table below).
##
##   demo("Ex_13b_rLMM_estimated_dispersion_known_vcov_BigWordClub_PartVI", package = "lmebayesCore")

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
##    groups yet). Mirrors demo("Ex_24_lmerb_dGamma_BigWordClub", package =
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
## 2. Design + priors: model_setup() / Prior_Setup_lmebayes() (first pass) /
##    per-group pwt_measurement calibration (2b) / Prior_Setup_lmebayes()
##    (final, tailored -- 2c) / pfamily_list()
##
## dispformula = ~school_id (matching the grouping factor exactly) requests
## Prior_Setup_lmebayes()'s per-group Block~1 calibration
## (ing_prior_measurement_group), consumed below by the Part VI extension
## (in place of Ex_13's direct dGamma_list() call).
##
## pwt_measurement is calibrated in two passes, same as Ex_13: 'ps_flat'
## below uses a flat scalar (0.1) purely as input to the per-group
## quantile-matching calibration in 2b/2c (fit_ref/Sigma_ranef/gamma_hat are
## all independent of pwt_measurement, so ps_flat's glmmTMB reference fit is
## identical to the final 'ps' produced in 2c). 'ps' -- the tailored,
## per-group pwt_measurement result of 2c -- is what the Part VI derivation
## (2d) and the sampler (Section 5) actually consume.
## ---------------------------------------------------------------------------
design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup (full-rank schools only) ===\n\n")
print(design)
stopifnot(all(design$re_rank))

## Same max_disp_perc = 0.8 as Ex_13 (Prior_Setup_lmebayes()'s own per-group
## sigma2_hat calibration is unchanged by Part VI -- only the window built
## from it, below, differs). pwt_measurement = 0.1 here is only the flat
## first-pass input to the per-group calibration in 2b/2c below -- it is not
## what the Part VI derivation/sampler end up using.
ps_flat <- Prior_Setup_lmebayes(
  form_lmer,
  data            = dat,
  pwt             = 0.01,
  dispformula     = ~school_id,
  max_disp_perc   = 0.8,
  pwt_measurement = 0.1
)
cat("\n=== Prior_Setup_lmebayes, first pass (flat pwt_measurement) ===\n\n")
print(ps_flat)

## group_name is not a formal on the routed export; attach it to 'group'
## instead of relying on substitute() (see .lmebayes_resolve_group_name()).
grp <- design$group
attr(grp, "group_name") <- design$group_name

group_levels <- levels(grp)
re_names     <- design$re_coef_names
p_re         <- length(re_names)

## ---------------------------------------------------------------------------
## 2b. Per-group pwt_measurement large enough to bring the pre-run
##     (noncentral, closed-form Gaussian plug-in) estimate of
##     pct_draws_outside below a target alpha_target -- see
##     inst/omega-ing-marginal-multivariate-t.md Section 9 and
##     data-raw/_scratch_group_pwt_measurement_noncentral.R (temporary
##     helper script, investigation-grade, not exported package code). Same
##     diagnostic as Ex_13 Section 2b -- see its longer note for the
##     rationale.
## ---------------------------------------------------------------------------
source("data-raw/_scratch_group_pwt_measurement_noncentral.R")
tab_pwt <- .tmp_group_pwt_measurement_noncentral(
  ps = ps_flat, design = design, group = grp, group_levels = group_levels,
  re_coef_names = re_names
)
.tmp_print_group_pwt_measurement_noncentral(tab_pwt)

w_pwt_vec <- stats::setNames(
  tab_pwt$w_star_floored[match(group_levels, tab_pwt$group)], group_levels
)
cat("\npwt_measurement vector (floored at 0.1), positional order matching levels(group):\n\n")
print(round(w_pwt_vec, 4))

## ---------------------------------------------------------------------------
## 2c. Feed w_pwt_vec into a SECOND, final Prior_Setup_lmebayes() call --
##     same formula/data/dispformula/max_disp_perc as 'ps_flat' above, but
##     with the tailored per-group pwt_measurement vector instead of the flat
##     0.1 scalar. This 'ps' -- not 'ps_flat' -- is what the Part VI
##     derivation (2d) and the sampler (Section 5) below consume.
## ---------------------------------------------------------------------------
ps <- Prior_Setup_lmebayes(
  form_lmer,
  data            = dat,
  pwt             = 0.01,
  dispformula     = ~school_id,
  max_disp_perc   = 0.8,
  pwt_measurement = w_pwt_vec
)
cat("\n=== Prior_Setup_lmebayes, final pass (tailored per-group pwt_measurement) ===\n\n")
print(ps)

rows_pwt_cmp <- character(0L)
for (lev in group_levels) {
  g1 <- ps_flat$ing_prior_measurement_group[[lev]]
  g2 <- ps$ing_prior_measurement_group[[lev]]
  rows_pwt_cmp <- c(rows_pwt_cmp, sprintf(
    "  %-6s  %8.4f  %8.4f  %12.4f  %12.4f  %10.4f  %10.4f\n",
    lev, w_pwt_vec[[lev]], g2$shape_ING, g1$rate, g2$rate, g1$disp_upper, g2$disp_upper
  ))
}
cat(
  "\n=== pwt_measurement = 0.1 (ps_flat) vs tailored (ps, final): shape/rate/disp_upper ===\n\n",
  sprintf("  %-6s  %8s  %8s  %12s  %12s  %10s  %10s\n",
          "group", "w_j", "shape_ING", "rate (flat)", "rate (final)",
          "disp_up(flat)", "disp_up(final)"),
  rows_pwt_cmp,
  sep = ""
)

## dNormal() Block~2 for every random-effect component: tau^2_k is *known*
## (fixed at its lmer REML estimate), so gamma_k has a conjugate Normal
## posterior -- no envelope/Gamma step. Also supplies Sigma_fixef below,
## the input the Part VI Omega_j extension propagates through W. Built from
## the TAILORED 'ps' (2c) -- pf itself is invariant to pwt_measurement, but
## this keeps it consistent with the ps that feeds everything else.
pf <- pfamily_list(ps)

## ---------------------------------------------------------------------------
## 2d. Part VI: per-group Block~1 measurement-dispersion prior with
##     Omega_j folded in (inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md Part VI).
##
## Built from the TAILORED 'ps' (2c) -- n_prior_j below (via
## ps$ing_prior_measurement_group[[lev]]$n_prior) picks up each group's own
## calibrated pwt_measurement automatically.
##
## Reproduces dGamma_list(ps, disp_upper_anchor = "symmetric")'s own two
## internal calls (.lmebayes_ing_prior_measurement_group_glm_inputs(),
## .lmebayes_compute_ing_prior_cal_from_sigma()) and its Part 0 Sigma_j
## shrinkage formula verbatim, adding only Sigma_j' = Sigma_j + Omega_j
## before calibration, and the same shape_w/rate_w symmetric-quantile
## window construction (Part II) applied to the resulting sigma2_hat,j'.
## ---------------------------------------------------------------------------
max_disp_perc <- 0.8
block_formula <- ps$block_formula
sd_tau        <- sqrt(diag(ps$Sigma_ranef))
re_names_all  <- design$re_coef_names
group_levels0 <- levels(dat$school_id)

part_vi_group <- stats::setNames(
  lapply(group_levels0, function(lev) {
    dat_j <- dat[dat$school_id == lev, , drop = FALSE]
    inp <- lmebayesCore:::.lmebayes_ing_prior_measurement_group_glm_inputs(
      lev = lev, dat_j = dat_j, block_formula = block_formula, sd_tau = sd_tau,
      family = gaussian(), intercept_source = "null_model", effects_source = "null_effects"
    )
    n_j          <- inp$n_j
    n_prior_j    <- ps$ing_prior_measurement_group[[lev]]$n_prior
    n_combined_j <- n_prior_j + n_j
    p_re         <- length(sd_tau)

    ## Part 0: coefficient-scale prior covariance from population sd_tau
    ## shrinkage (verbatim from .lmebayes_calibrate_ing_prior_measurement_
    ## group(), R/mixed_rmerb_helpers.R).
    pwt_j <- diag(inp$V0)
    pwt_j <- pwt_j / (pwt_j + inp$sd_vec^2)
    if (length(pwt_j) == 1L) {
      Sigma_j <- ((1 - pwt_j) / pwt_j) * inp$V0
    } else {
      scale_vec <- sqrt((1 - pwt_j) / pwt_j)
      Sigma_j <- inp$V0 * outer(scale_vec, scale_vec)
    }

    ## Part VI: model-derived Omega_j, diagonal across RE components (each
    ## gamma_k calibrated independently in pf/ps$prior_list).
    Omega_j <- matrix(0, nrow = length(inp$var_names), ncol = length(inp$var_names),
                       dimnames = list(inp$var_names, inp$var_names))
    for (k in re_names_all) {
      Wk_row <- design$W[[k]][lev, , drop = FALSE]
      Sigma_fixef_k <- ps$prior_list[[k]]$Sigma_fixef
      Omega_j[k, k] <- as.numeric(Wk_row %*% Sigma_fixef_k %*% t(Wk_row))
    }

    cal <- lmebayesCore:::.lmebayes_compute_ing_prior_cal_from_sigma(
      inp, Sigma_j + Omega_j, n_prior_j
    )

    shape_w <- (n_combined_j + 1) / 2 + p_re / 2
    rate_w  <- cal$dispersion * (n_combined_j + p_re - 1) / 2
    disp_lower <- 1 / qgamma(max_disp_perc,     shape = shape_w, rate = rate_w)
    disp_upper <- 1 / qgamma(1 - max_disp_perc, shape = shape_w, rate = rate_w)

    list(
      sigma2_hat = cal$dispersion,
      shape      = cal$shape_ING,
      rate       = cal$rate,
      disp_lower = disp_lower,
      disp_upper = disp_upper,
      omega_j    = Omega_j
    )
  }),
  group_levels0
)

cat("\n=== Part VI per-group sigma^2_j: model-derived Omega_j vs Ex_13's symmetric-only prior ===\n\n")
part_vi_tab <- do.call(rbind, lapply(group_levels0, function(lev) {
  g <- part_vi_group[[lev]]
  data.frame(
    group      = lev,
    sigma2_hat = g$sigma2_hat,
    disp_lower = g$disp_lower,
    disp_upper = g$disp_upper
  )
}))
num_cols <- vapply(part_vi_tab, is.numeric, logical(1L))
part_vi_tab[num_cols] <- lapply(part_vi_tab[num_cols], round, digits = 3)
print(part_vi_tab, row.names = FALSE)

## ---------------------------------------------------------------------------
## 3. Arguments matrix_args_lmm() would build for rlmerb() -- assembled here
##    by hand so rLMMindepNormalGamma_reg_known_vcov() can be called directly.
##
## The routed export's 'prior_list' for a per-group ING Block~1 dispersion is
## NOT the dGamma() pfamily list itself -- it is a flat list with 'mu'/'Sigma'
## (the Block~2 hyperparameter prior, same shape .rLMM_validate_ing_
## measurement_prior_list() expects for the fixed-vcov/estimated-vcov cases)
## plus 'shape_group'/'rate_group'/'disp_lower_group'/'disp_upper_group'
## (one named-by-group-level numeric vector each) -- extracted here from
## part_vi_group (Section 2d) instead of a dGamma() pfamily list, mirroring
## .lmebayes_resolve_dispersion_ranef_group_list() /
## .lmebayes_ing_measurement_prior_list_group() in mixed_rmerb_helpers.R.
## ---------------------------------------------------------------------------

shape_group      <- stats::setNames(numeric(length(group_levels)), group_levels)
rate_group       <- stats::setNames(numeric(length(group_levels)), group_levels)
disp_lower_group <- stats::setNames(numeric(length(group_levels)), group_levels)
disp_upper_group <- stats::setNames(numeric(length(group_levels)), group_levels)
for (lev in group_levels) {
  g <- part_vi_group[[lev]]
  shape_group[[lev]]      <- g$shape
  rate_group[[lev]]       <- g$rate
  disp_lower_group[[lev]] <- g$disp_lower
  disp_upper_group[[lev]] <- g$disp_upper
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
  "\n=== Per-group sigma^2_j ING prior (Part VI): %d groups, sigma^2_hat range [%.4f, %.4f] ===\n\n",
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
## 5. Direct call: rLMMindepNormalGamma_reg_known_vcov()
##
## Two-block Gibbs with an ING Block~1 sweep (sigma^2_j is estimated, so the
## joint posterior is not exactly Gaussian even though Block~2 is dNormal()):
## an optional pilot stage (gap_tol/mode_gap_max) recenters the main stage's
## starting point away from the ICM mode, then Theorem~3 calibrates the
## number of inner sweeps per stored draw.
##
## progbar/verbose match demo("Ex_24_lmerb_dGamma_BigWordClub", package =
## "lmebayes"): that demo calls lmerb() without overriding progbar/verbose,
## and lmerb()'s own formals are progbar = NULL (falsy -- no bar shown) and
## a hardcoded verbose = TRUE passed to rlmerb().
## ---------------------------------------------------------------------------
fit <- rLMMindepNormalGamma_reg_known_vcov(
  n            = 1000L,
  y            = design$y,
  D            = design$D,
  group        = grp,
  W            = design$W,
  prior_list   = prior_list,
  pfamily_list = pf,
  progbar      = FALSE,
  verbose      = TRUE
)

source("data-raw/_scratch_rss_ellipsoid_test.R", local = FALSE)  # defines .tmp_rss_ellipsoid_test only if you comment out/skip the run_one() calls at the bottom

tab_13b <- .tmp_rss_ellipsoid_test(
  fit           = fit,
  D             = design$D,
  y             = design$y,
  group         = grp,
  group_name    = design$group_name,
  re_coef_names = re_names,
  shape_group   = shape_group,
  rate_group    = rate_group
)
print(tab_13b[order(tab_13b$p_value), ], row.names = FALSE, digits = 4)

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
## 6. sigma^2_j: post-sweep draws' per-group means vs Part VI-calibrated
##    prior mean, and vs the pooled lmer REML sigma^2.
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
  "\n=== sigma^2_j: post mean vs Part VI-calibrated prior mean (all groups) ===\n\n",
  sprintf("  %-6s  %10s  %10s  %10s\n",
          "group", "post mean", "prior mean", "[window]"),
  rows_disp,
  sep = ""
)

## ---------------------------------------------------------------------------
## 7. Block 2 fixed effects: Gibbs (MCMC) vs ICM mean vs glmmTMB fixef
##    (+ uncertainty)
##
## No exact-iid engine exists for this model (sigma^2_j is estimated per
## group, so the joint posterior is not exactly Gaussian) -- 'gibbs mean'/
## 'gibbs SD' (the lmebayesCore output: posterior mean/SD of gamma_k across
## the main-stage MCMC draws) is compared directly to the same
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

## ---------------------------------------------------------------------------
## 7b. Extended (Omega)-aware local rate diagnostic
##     (inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 14-15)
##
## Same diagnostic as Ex_13 Section 7b: omega_ing$omega/$e are evaluated at
## the now-*completed* sampler's own posterior-mean state --
## 1/fit$dispersion_ranef.mean[[lev]] for Omega_j, and e_j = y_j - D_j %*%
## (posterior-mean beta_j from fit$coefficients) -- rather than the
## certified corner's disp_upper_group/fit_ref (glmmTMB) combination, so
## both quantities come from the SAME fitted object and describe one
## mutually-consistent reference state (see Ex_13 Section 7b's longer note
## for why mixing corner/external-reference-fit quantities can put Omega_j
## and e_j out of step for a poorly-predicted group). This is a
## *diagnostic*, not a certified bound: see two_block_rate_ing()'s own
## documentation for why no analogue of Theorem~3/Corollary~1 exists for
## the extended chain.
## ---------------------------------------------------------------------------
n_group  <- stats::setNames(as.numeric(table(grp)), group_levels)
e_post   <- lmebayesCore:::.lmebayes_posterior_group_residuals(
  fit, y = design$y, D = design$D, group = grp,
  group_name = design$group_name, re_coef_names = re_names
)
omega_post <- stats::setNames(
  1 / fit$dispersion_ranef.mean[group_levels], group_levels
)

omega_ing <- list(
  scope = "per_group",
  omega = omega_post,
  shape = shape_group[group_levels],
  n     = n_group,
  e     = e_post
)

prior_list_block1_rate <- list(
  Sigma      = as.matrix(ps$Sigma_ranef),
  dispersion = disp_upper_group[group_levels]
)
prior_list_block2_rate <- lapply(pf, function(pfk) {
  pl <- pfk$prior_list
  list(mu = pl$mu, Sigma = pl$Sigma, dispersion = pl$dispersion)
})

rate_ext <- two_block_rate_ing(
  x = design$D, block = grp, x_hyper = design$W,
  prior_list_block1 = prior_list_block1_rate,
  prior_list_block2 = prior_list_block2_rate,
  omega_ing = omega_ing
)
cat("\n=== Extended (Omega)-aware local rate diagnostic (Section 14-15) ===\n\n")
print(rate_ext)

## ---------------------------------------------------------------------------
## 7c. Empirical worst-case: same (Omega)-aware extended rate, evaluated at
##     EVERY main-stage draw instead of the single posterior-mean reference
##     state above. Same diagnostic as Ex_13 Section 7c -- see its longer
##     note for the rationale (pilot-draw pmax() scan analogue, applied to
##     the extended system and the main-stage output; still a diagnostic,
##     not a certified bound).
## ---------------------------------------------------------------------------
omega_spec <- list(
  scope = "per_group",
  shape = shape_group[group_levels],
  n     = n_group
)

rate_emp <- lmebayesCore:::.two_block_rate_ing_over_draws(
  fit = fit, n_draws = n_draws,
  x = design$D, block = grp, x_hyper = design$W,
  prior_list_block1 = prior_list_block1_rate,
  prior_list_block2 = prior_list_block2_rate,
  group_name = design$group_name, re_coef_names = re_names,
  y = design$y, D = design$D,
  omega_spec = omega_spec
)

cat(sprintf(
  paste0(
    "\n=== Empirical (Omega)-aware local rate over n = %d main-stage draws ===\n\n",
    "  lambda* (extended, single reference state)     = %.6f\n",
    "  lambda* (extended, empirical max over %d draws) = %.6f  (draw #%d)\n",
    "  draws with lambda*(extended) >= 1: %d / %d\n"
  ),
  rate_emp$n_draws,
  rate_ext$lambda_star,
  rate_emp$n_draws, rate_emp$lambda_star_max, rate_emp$i_max,
  rate_emp$n_over_one, rate_emp$n_draws
))

## ---------------------------------------------------------------------------
## 7d. SCRATCH: Omega_j-MARGINALIZED (Section 16.2) lambda_star, also
##     evaluated at every main-stage draw -- see
##     data-raw/_scratch_lambda_star_marginal_over_draws.R (temporary,
##     investigation only, not package code). Same rationale as Ex_13
##     Section 7d: unlike 7c above (which plugs each draw's
##     independently-sampled (beta_j, Omega_j) pair into the *joint* Section
##     14 Hessian), this integrates Omega_j out analytically first (Section
##     16), so each draw's effective measurement precision is always the one
##     implied by that SAME draw's own beta_j. Draws where some group's
##     Lambda + H_j(beta_j) block isn't PD are skipped (flagged per-group)
##     rather than forced through.
## ---------------------------------------------------------------------------
source("data-raw/_scratch_lambda_star_marginal_over_draws.R")
inp_marg <- lmebayesCore:::.two_block_rate_inputs(
  x = design$D, block = grp, x_hyper = design$W,
  prior_list_block1 = prior_list_block1_rate,
  prior_list_block2 = prior_list_block2_rate
)
blocks_marg <- lmebayesCore:::.two_block_S_P11(inp_marg)
group_setup_marg <- .tmp_marginal_group_setup(
  D = design$D, y = design$y, group = grp, group_levels = group_levels,
  re_coef_names = re_names, shape_group = shape_group, rate_group = rate_group
)
res_marg <- .tmp_lambda_star_marginal_over_draws(
  fit = fit, n_draws = n_draws, y = design$y,
  group_name = design$group_name, group_setup = group_setup_marg,
  inp = inp_marg, blocks = blocks_marg
)
.tmp_print_marginal_over_draws_summary(res_marg)

## ---------------------------------------------------------------------------
## 7e. SCRATCH: split-support revised end-of-simulation TV bound
##     (inst/BLOCK_GIBBS_ERGODICITY_ING.md Section 17) -- see
##     data-raw/_scratch_tv_bound_revised.R (temporary, investigation only,
##     not package code). Reuses res_marg (7d) directly: A^c is the set of
##     main-stage draws where some group's Lambda + H_j wasn't PD, or the
##     draw's own lambda_star_marginal exceeded lambda_star_used -- the
##     lambda_star that actually calibrated this run's m_convergence (the
##     marginal safeguard if it was valid for this pilot, else the base/
##     extended upper bound it fell back to).
## ---------------------------------------------------------------------------
source("data-raw/_scratch_tv_bound_revised.R")
lambda_star_used <- if (isTRUE(fit$convergence_info$marginal_rate_valid)) {
  fit$convergence_info$lambda_star_marginal
} else {
  fit$convergence_info$lambda_star_upper
}
tv_revised <- .tmp_tv_bound_revised(
  res_marg, lambda_star_used, tv_tol = fit$convergence_info$tv_tol
)
.tmp_print_tv_bound_revised(tv_revised)

## Combined mean-bias/Var_final-ratio charts (Claims 1 and 3 of the two-block
## Gibbs ergodicity reference): rLMMindepNormalGamma_reg_known_vcov() goes
## through the sweeps-outer/chains-inner pilot/main engine, so both
## fit$pilot$sweep_history and fit$sweep_history carry cov_by_sweep (as does
## rLMMNormal_reg_known_vcov(sim_method = "TWO_BLOCK_GIBBS")'s single main
## stage now that its engine is rGLMM_sweep()-based too). Per-group dispersion
## is *estimated* here (that is the point of this demo), so
## plot_mean_convergence()/plot_var_convergence()'s fit-object methods
## automatically find no exact reference available and fall back to the
## empirical mean_final/Var_final (last-sweep cross-chain mean/covariance).
## 'stage' selects fit$pilot$sweep_history ("pilot") or fit$sweep_history
## ("main"); n_chains defaults to the correct chain count for each stage
## (fit$pilot_chisq$n_pilot vs fit$n). The pilot stage's whole job is to
## *recenter* the main stage's starting point away from the ICM mode, so the
## mean-bias chart is the more direct diagnostic there; the main stage is
## where the variance/uncertainty calibration check matters most -- both
## charts are shown for both stages for consistency.
for (stg in c("pilot", "main")) {
  plot_mean_convergence(fit, whitened = FALSE, stage = stg)
  plot_mean_convergence(fit, whitened = TRUE, stage = stg)
  plot_var_convergence(fit, whitened = FALSE, stage = stg)
  plot_var_convergence(fit, whitened = TRUE, stage = stg)
}

## ---------------------------------------------------------------------------
## 8. Sweep-history diagnostics: cross-chain mean/SD vs inner sweep, for
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
## 9. Random effects: MCMC mean (per group, per draw average) vs ICM mode
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
