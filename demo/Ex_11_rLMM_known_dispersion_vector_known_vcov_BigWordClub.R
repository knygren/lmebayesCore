## Demo: rLMMNormal_reg_known_vcov() called directly on bayesrules::big_word_club
## with a per-group KNOWN (fixed, not estimated) observation dispersion vector
##
## Case 2 of 5: KNOWN observation dispersion, but *per-group* -- sigma^2_j
## varies by school_id and is fixed/given for the duration of sampling (NOT
## a dGamma() prior to be sampled -- contrast with
## demo("Ex_13_rLMM_estimated_dispersion_known_vcov_BigWordClub", package =
## "lmebayesCore"), where sigma^2_j is ESTIMATED per group). Also contrast
## with demo("Ex_10_rLMM_known_dispersion_known_vcov_BigWordClub", package =
## "lmebayesCore"), where sigma^2 is a single value pooled across all
## groups. KNOWN random-effect variance components as in Ex_10 (every
## Block~2 pfamily component is dNormal(), so tau^2_k is fixed, not
## sampled).
##
## rLMMNormal_reg_known_vcov()'s 'prior_list$dispersion' accepts either a
## single positive scalar (Ex_10) or a length-J named vector, one known
## value per group level (.rLMM_validate_fixed_dispersion_vector();
## equivalently lmerb()'s 4th 'dispersion_ranef' mode, "fixed_vector", see
## .lmebayes_resolve_dispersion_ranef_fixed_vector() in
## mixed_rmerb_helpers.R) -- this demo exercises that vector form.
##
## The per-group sigma^2_j values used here are
## Prior_Setup_lmebayes(..., dispformula = ~school_id)'s 'sigma2_group'
## field: per-group point estimates read off a glmmTMB heteroscedastic
## reference fit's dispersion linear predictor. Prior_Setup_lmebayes()
## documents 'sigma2_group' as "diagnostic only -- not the value fed to the
## sampler" (its own sampler-facing use is the *estimated* per-group
## dGamma_list() route, Ex_13/Ex_14) -- this demo deliberately repurposes it
## as a *known, fixed* input instead, and (per the request that motivated
## this demo) runs both rLMMNormal_reg_known_vcov() engines
## (sim_method = "DEFAULT", the exact-iid route, and sim_method =
## "TWO_BLOCK_GIBBS") and compares both against the same glmmTMB reference
## fit's fixed and random effects.
##
##   demo("Ex_11_rLMM_known_dispersion_vector_known_vcov_BigWordClub", package = "lmebayesCore")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}
if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop("This demo requires the 'glmmTMB' package.", call. = FALSE)
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
## 1. Prior_Setup_lmebayes(dispformula = ~school_id)'s per-group calibration
##    (ing_prior_measurement_group, used for sigma2_group here) runs a
##    within-group regression for every group and requires each group's Z_j
##    to be full column rank (unlike the *engine*'s fixed_vector dispersion
##    mode itself, which has no such requirement -- see
##    .lmebayes_resolve_dispersion_ranef_fixed_vector()'s doc comment). No
##    ING accept/reject envelope is built for this route, so -- unlike
##    Ex_13/Ex_14 -- school_id 2/18 do not need to be excluded here.
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

## ---------------------------------------------------------------------------
## 2. Design + priors: model_setup() / Prior_Setup_lmebayes()
##
## dispformula = ~school_id (matching the grouping factor exactly) requests
## Prior_Setup_lmebayes()'s glmmTMB-based per-group calibration, which is
## also where Block~2's prior_list/Sigma_ranef/tau^2_k plug-ins come from
## here (calibration_source == "glmmTMB").
## ---------------------------------------------------------------------------
design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup (full-rank schools only) ===\n\n")
print(design)
stopifnot(all(design$re_rank))

ps <- Prior_Setup_lmebayes(
  form_lmer,
  data        = dat,
  pwt         = 0.01,
  dispformula = ~school_id
)
cat("\n=== Prior_Setup_lmebayes (per-group Block~1 calibration) ===\n\n")
print(ps)

## dNormal() Block~2 for every random-effect component: tau^2_k is *known*
## (fixed at its glmmTMB REML plug-in), so gamma_k has a conjugate Normal
## posterior -- no envelope/Gamma step.
pf <- pfamily_list(ps)

## ---------------------------------------------------------------------------
## 3. Arguments assembled by hand for rLMMNormal_reg_known_vcov(): 'group' with
##    an attached 'group_name' (matrix_args_lmm()'s usual pattern), and a
##    per-group KNOWN dispersion vector.
## ---------------------------------------------------------------------------
grp <- design$groups
attr(grp, "group_name") <- design$group_name
group_levels <- levels(grp)
re_names     <- design$re_coef_names

## sigma2_group is a diagnostic-only field in Prior_Setup_lmebayes(); here we
## deliberately treat it as KNOWN/fixed and pass it straight through as
## prior_list$dispersion's per-group vector (not sampled at all, unlike
## Ex_13/Ex_14's dGamma_list() route).
disp_known <- ps$sigma2_group[group_levels]
prior_list <- list(dispersion = disp_known)

cat(sprintf(
  "\n=== Known per-group sigma^2_j (from glmmTMB dispformula fit): range [%.4f, %.4f] ===\n\n",
  min(disp_known), max(disp_known)
))

## ---------------------------------------------------------------------------
## 4. glmmTMB reference fit (the *same* fit Prior_Setup_lmebayes() calibrated
##    Block~2/sigma2_group from -- no re-fitting needed).
## ---------------------------------------------------------------------------
fit_ref <- ps$fit_ref
cat("\n=== glmmTMB reference fit (dispformula = ~school_id) ===\n\n")
print(summary(fit_ref))

## NOTE on the "Residual  NA  NA" row in the VarCorr table above: with
## dispformula = ~school_id, the observation-level residual variance is no
## longer a single scalar -- it is the output of glmmTMB's own dispersion
## GLM (log link by default for family = gaussian()), so there is no single
## number to put in that summary slot and glmmTMB prints NA. The actual
## per-group values ARE in this same summary, in two equivalent forms:
##  (a) the "Dispersion model:" coefficient table further above gives each
##      school's *log*-dispersion as a contrast from the reference school
##      (first factor level) -- e.g. its (Intercept) row is
##      log(sigma^2) for the reference school, and each school_id<k> row is
##      that school's *offset* from it on the log scale;
##  (b) Section 8 below prints the same values already exponentiated back
##      to the natural (sigma^2) scale for every school -- this is exactly
##      predict(fit_ref, type = "disp"), averaged within each group, which
##      is also precisely what ps$sigma2_group holds and what 'disp_known'
##      (prior_list$dispersion) above was built from.
## As a spot check: Dispersion model (Intercept) = 3.37823 on the log scale
## for the reference school ('2') exponentiates to exp(3.37823) = 29.32,
## matching sigma^2_2 in Section 8 below.

## ---------------------------------------------------------------------------
## 5. Direct calls: rLMMNormal_reg_known_vcov() with BOTH engines --
##    sim_method = "DEFAULT" (exact iid draws from the closed-form joint
##    Gaussian posterior) and sim_method = "TWO_BLOCK_GIBBS" (Theorem~3-
##    calibrated Gibbs sweeps) -- both target the exact same posterior here
##    (fully Gaussian-conjugate: known per-group dispersion, known vcov).
##
## progbar/verbose match demo("Ex_12_lmerb_BigWordClub", package =
## "lmebayes"): that demo calls lmerb() without overriding progbar/verbose,
## and lmerb()'s own formals are progbar = NULL (falsy -- no bar shown) and
## a hardcoded verbose = TRUE passed to rlmerb().
## ---------------------------------------------------------------------------
fit_iid <- rLMMNormal_reg_known_vcov(
  n            = 10000L,
  y            = design$y,
  D            = design$Z,
  group        = grp,
  W            = design$X_hyper,
  prior_list   = prior_list,
  pfamily_list = pf,
  progbar      = FALSE,
  verbose      = TRUE,
  sim_method   = "DEFAULT"
)
cat(sprintf("\nsim_method_used (fit_iid): %s\n", fit_iid$sim_method_used))

fit_gibbs <- rLMMNormal_reg_known_vcov(
  n            = 10000L,
  y            = design$y,
  D            = design$Z,
  group        = grp,
  W            = design$X_hyper,
  prior_list   = prior_list,
  pfamily_list = pf,
  progbar      = FALSE,
  verbose      = TRUE,
  sim_method   = "TWO_BLOCK_GIBBS"
)
cat(sprintf("sim_method_used (fit_gibbs): %s (m_convergence = %d)\n",
            fit_gibbs$sim_method_used, fit_gibbs$m_convergence))

n_draws <- nrow(fit_iid$fixef[[re_names[1L]]])

## ---------------------------------------------------------------------------
## 6. Block 2 fixed effects: DEFAULT (iid) vs TWO_BLOCK_GIBBS draws means vs
##    exact ICM mean vs glmmTMB fixef, alongside BOTH sides' own uncertainty
##    ('iid SD' = posterior SD of gamma_k from the draws, i.e. the Bayesian
##    analogue of glmmTMB's 'Std.Error'; NOT the Monte Carlo SE(mean), which
##    is much smaller and only measures how precisely n = 1000 draws pin down
##    the posterior mean itself).
##
## dispformula = ~school_id spends 34 parameters (1 intercept + 33 school_id
## deviations) estimating per-group dispersion from only 400 observations,
## leaving several fixed effects only weakly identified by glmmTMB itself
## (e.g. private_school: z = -0.11; distracted_a1: z = 0.41) -- comparing
## raw point estimates alone (as an earlier draft of this demo did) can look
## alarming purely from that weak identification, not from any actual
## disagreement; the 'diff (SE units)' column re-expresses the iid-mean vs
## glmmTMB gap in units of glmmTMB's own Std.Error, which is the right scale
## to judge it on.
## ---------------------------------------------------------------------------
fe_ref  <- lmebayesCore:::.lmebayes_reference_fixef(fit_ref)
se_ref  <- sqrt(diag(lmebayesCore:::.lmebayes_reference_vcov(fit_ref)))

## Build every row first (no printing) so the header/dashes/rows below print
## as one contiguous block -- otherwise demo()'s echo = TRUE interleaves this
## loop's own source between the header and the first data row.
rows_fe <- character(0L)
for (k in re_names) {
  dm_iid   <- colMeans(fit_iid$fixef[[k]])
  sd_iid   <- apply(fit_iid$fixef[[k]], 2L, sd)
  dm_gibbs <- colMeans(fit_gibbs$fixef[[k]])
  sd_gibbs <- apply(fit_gibbs$fixef[[k]], 2L, sd)
  icm_k    <- fit_iid$fixef.mode[[k]]
  for (nm in names(dm_iid)) {
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
    diff_se <- (dm_iid[[nm]] - fe_val) / se_val
    rows_fe <- c(rows_fe, sprintf(
      "  %-18s  %-28s  %10.4f  %8.4f  %10.4f  %8.4f  %10.4f  %10.4f  %8.4f  %8.2f\n",
      k, nm, dm_iid[[nm]], sd_iid[[nm]], dm_gibbs[[nm]], sd_gibbs[[nm]],
      icm_k[[nm]], fe_val, se_val, diff_se
    ))
  }
}

## Single cat() call for the whole block (title + header + dashes + rows +
## footnote): demo()'s echo = TRUE prints one "> ..." per top-level R
## statement, so splitting the table across multiple cat() calls -- even
## adjacent ones -- still interleaves an echoed statement between each
## piece. One call in, one block out.
cat(
  "\n=== Block 2 fixed effects: iid vs Gibbs vs ICM mean vs glmmTMB fixef (+ uncertainty) ===\n\n",
  sprintf("  %-18s  %-28s  %10s  %8s  %10s  %8s  %10s  %10s  %8s  %8s\n",
          "RE component", "parameter", "iid mean", "iid SD", "gibbs mean", "gibbs SD",
          "ICM mean", "glmmTMB", "glmm SE", "diff(SE)"),
  sprintf("  %-18s  %-28s  %10s  %8s  %10s  %8s  %10s  %10s  %8s  %8s\n",
          strrep("-", 18L), strrep("-", 28L), strrep("-", 10L), strrep("-", 8L),
          strrep("-", 10L), strrep("-", 8L), strrep("-", 10L), strrep("-", 10L),
          strrep("-", 8L), strrep("-", 8L)),
  rows_fe,
  "\n  diff(SE) = (iid mean - glmmTMB estimate) / glmmTMB Std.Error -- |diff(SE)| < ~1-2\n",
  "  is well within glmmTMB's own uncertainty for that coefficient, not a discrepancy.\n",
  sep = ""
)

## ---------------------------------------------------------------------------
## 7. Random effects: glmmTMB "full" coefficient vs sampler ranef.mode (both
##    engines), all groups -- same column layout as Section 6's fixed-effects
##    table (iid mean/SD, gibbs mean/SD, ICM mean, glmmTMB, diff(SE)), except
##    glmmTMB has no per-group random-effect standard error to report here
##    (unlike the fixed-effects case, where fixef()'s Std.Error comes
##    straight out of vcov(fit_ref)) -- 'iid SD' (the posterior SD of beta_j
##    across the iid engine's own draws) is used as the uncertainty scale for
##    'diff(SE)' instead.
##
## Mirrors Ex_10's Section 7 (build_mu_all()-based reconstruction, needed
## because glmmTMB::coef()[[k]], like lme4::coef()[[k]], omits any cross-
## level covariate's contribution to RE component k whenever W_k has
## non-intercept columns), but sourced from the glmmTMB reference fit_ref
## instead of an lme4::lmer() fit.
## ---------------------------------------------------------------------------
grp_col  <- design$group_name
grp_levs <- rownames(fit_iid$ranef.mode)

## Per-group, per-RE-component mean/SD of beta_j across each engine's own
## draws (fit_*$coefficients: long data.frame, one row per (draw, group)) --
## same construction as Ex_10/Ex_12/Ex_13/Ex_14's "MCMC mean vs ICM" tables.
re_draws_mean_iid <- tapply(
  seq_len(nrow(fit_iid$coefficients)),
  fit_iid$coefficients[[grp_col]],
  function(idx) colMeans(fit_iid$coefficients[idx, re_names, drop = FALSE]),
  simplify = FALSE
)
re_draws_sd_iid <- tapply(
  seq_len(nrow(fit_iid$coefficients)),
  fit_iid$coefficients[[grp_col]],
  function(idx) apply(fit_iid$coefficients[idx, re_names, drop = FALSE], 2L, sd),
  simplify = FALSE
)
re_draws_mean_gibbs <- tapply(
  seq_len(nrow(fit_gibbs$coefficients)),
  fit_gibbs$coefficients[[grp_col]],
  function(idx) colMeans(fit_gibbs$coefficients[idx, re_names, drop = FALSE]),
  simplify = FALSE
)
re_draws_sd_gibbs <- tapply(
  seq_len(nrow(fit_gibbs$coefficients)),
  fit_gibbs$coefficients[[grp_col]],
  function(idx) apply(fit_gibbs$coefficients[idx, re_names, drop = FALSE], 2L, sd),
  simplify = FALSE
)

coef_raw <- as.data.frame(lmebayesCore:::.lmebayes_reference_coef(fit_ref)[[grp_col]])
if (!identical(names(coef_raw), re_names) && all(re_names %in% names(coef_raw))) {
  coef_raw <- coef_raw[, re_names, drop = FALSE]
}

.fe_name_for_ref <- function(k, col, fe) {
  if (k == "(Intercept)") {
    if (col %in% names(fe)) col else NA_character_
  } else if (col == "(Intercept)") {
    if (k %in% names(fe)) k else NA_character_
  } else {
    cand <- c(paste0(col, ":", k), paste0(k, ":", col))
    hit  <- cand[cand %in% names(fe)]
    if (length(hit)) hit[1L] else NA_character_
  }
}

fixef_ref <- lapply(re_names, function(k) {
  cols_k <- colnames(design$X_hyper[[k]])
  fe_nms <- vapply(cols_k, .fe_name_for_ref, character(1L), k = k, fe = fe_ref)
  miss <- is.na(fe_nms) | !fe_nms %in% names(fe_ref)
  if (any(miss)) {
    stop(
      "glmmTMB fixef missing term(s) for W[[", k, "]]: ",
      paste(cols_k[miss], collapse = ", "),
      call. = FALSE
    )
  }
  mu_k <- vapply(fe_nms, function(nm) unname(fe_ref[nm]), numeric(1L))
  names(mu_k) <- cols_k
  mu_k
})
names(fixef_ref) <- re_names

coef_anchor <- vapply(re_names, function(k) {
  if (k == "(Intercept)") unname(fe_ref["(Intercept)"]) else unname(fe_ref[k])
}, numeric(1L))

mu_all_ref <- build_mu_all(design, fixef_ref)$mu_all

## Build every row first (no printing) so the header/dashes/rows below print
## as one contiguous block -- otherwise demo()'s echo = TRUE interleaves this
## loop's own source between the header and the first data row.
rows_re <- character(0L)
for (lev in grp_levs) {
  glmm_k <- vapply(re_names, function(k) {
    mu_all_ref[k, lev] + (coef_raw[lev, k] - coef_anchor[[k]])
  }, numeric(1L))
  icm_k <- unname(fit_iid$ranef.mode[lev, re_names])
  for (i in seq_along(re_names)) {
    k        <- re_names[[i]]
    dm_iid   <- re_draws_mean_iid[[lev]][[k]]
    sd_iid   <- re_draws_sd_iid[[lev]][[k]]
    dm_gibbs <- re_draws_mean_gibbs[[lev]][[k]]
    sd_gibbs <- re_draws_sd_gibbs[[lev]][[k]]
    diff_se  <- (dm_iid - glmm_k[[i]]) / sd_iid
    rows_re <- c(rows_re, sprintf(
      "  %-6s  %-18s  %10.4f  %8.4f  %10.4f  %8.4f  %10.4f  %10.4f  %8.2f\n",
      lev, k, dm_iid, sd_iid, dm_gibbs, sd_gibbs, icm_k[[i]], glmm_k[[i]], diff_se
    ))
  }
}

## Single cat() call for the whole block (title + header + dashes + rows +
## footnote): demo()'s echo = TRUE prints one "> ..." per top-level R
## statement, so splitting the table across multiple cat() calls -- even
## adjacent ones -- still interleaves an echoed statement between each
## piece. One call in, one block out.
cat(
  "\n=== Random effects: iid vs Gibbs vs ICM mean vs glmmTMB 'full' coefficient (+ uncertainty) ===\n\n",
  sprintf("  %-6s  %-18s  %10s  %8s  %10s  %8s  %10s  %10s  %8s\n",
          "group", "RE component", "iid mean", "iid SD", "gibbs mean", "gibbs SD",
          "ICM mean", "glmmTMB", "diff(SE)"),
  sprintf("  %-6s  %-18s  %10s  %8s  %10s  %8s  %10s  %10s  %8s\n",
          strrep("-", 6L), strrep("-", 18L), strrep("-", 10L), strrep("-", 8L),
          strrep("-", 10L), strrep("-", 8L), strrep("-", 10L), strrep("-", 10L),
          strrep("-", 8L)),
  rows_re,
  "\n  glmmTMB = build_mu_all(design, fixef_ref)$mu_all + (coef(fit_ref) - coef_anchor).\n",
  "  diff(SE) = (iid mean - glmmTMB) / iid SD -- glmmTMB has no published per-group\n",
  "  random-effect standard error to compare against (unlike Section 6's fixed\n",
  "  effects, whose 'glmm SE' comes straight from vcov(fit_ref)), so the iid\n",
  "  engine's own posterior SD of beta_j is used as the uncertainty scale instead.\n",
  sep = ""
)

## ---------------------------------------------------------------------------
## 8. Known per-group sigma^2_j: the exact input fed to both samplers,
##    reproduced here for reference (identical to Prior_Setup_lmebayes()'s
##    printed 'ps$sigma2_group', by construction -- shown for completeness,
##    not a new estimate).
## ---------------------------------------------------------------------------
cat("\n=== Known per-group sigma^2_j fed to both samplers (all groups) ===\n\n")
for (lev in group_levels) {
  cat(sprintf("  %-6s  sigma^2_%s = %.4f\n", lev, lev, disp_known[[lev]]))
}
