## Demo: rLMMNormal_reg_known_vcov() called directly on bayesrules::big_word_club
##
## Case 1 of 5: KNOWN, single pooled observation dispersion (sigma^2 fixed
## at its lmer REML estimate, shared across all groups) and KNOWN
## random-effect variance components (every Block~2 pfamily component is
## dNormal(), so tau^2_k is fixed, not sampled). Compare
## demo("Ex_11_rLMM_known_dispersion_vector_known_vcov_BigWordClub", package
## = "lmebayesCore") for the same vcov case with a *per-group* known
## (fixed, not estimated) dispersion vector instead of a pooled scalar.
##
## Same model as demo("Ex_12_lmerb_BigWordClub", package = "lmebayes"), but
## this script calls rLMMNormal_reg_known_vcov() directly instead of going
## through lmerb()/rlmerb(): model_setup(), Prior_Setup_lmebayes(), and
## pfamily_list() (all exported from lmebayesCore) build the design and
## priors, then the script assembles by hand the exact 'group'/'prior_list'
## arguments that matrix_args_lmm() builds internally for rlmerb(), and calls
## the matrix-level export directly.
##
##   demo("Ex_10_rLMM_known_dispersion_known_vcov_BigWordClub", package = "lmebayesCore")

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
## 1. Design + priors: model_setup() / Prior_Setup_lmebayes() / pfamily_list()
## ---------------------------------------------------------------------------
design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup ===\n\n")
print(design)

ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)
cat("\n=== Prior_Setup_lmebayes ===\n\n")
print(ps)

## dNormal() Block~2 for every random-effect component: tau^2_k is *known*
## (fixed at its lmer REML estimate), so gamma_k has a conjugate Normal
## posterior -- no envelope/Gamma step.
pf <- pfamily_list(ps)

## ---------------------------------------------------------------------------
## 2. Arguments matrix_args_lmm() would build for rlmerb() -- assembled here
##    by hand so rLMMNormal_reg_known_vcov() can be called directly.
## ---------------------------------------------------------------------------

## The routed export has no 'group_name' formal; attach it to 'group' instead
## of relying on substitute() (see .lmebayes_resolve_group_name()).
grp <- design$groups
attr(grp, "group_name") <- design$group_name

## Known observation dispersion: fixed sigma^2 from the lmer REML fit.
prior_list <- list(dispersion = ps$dispersion_ranef)

cat(sprintf(
  "\n=== Known observation dispersion: sigma^2 = %.4f (lmer REML) ===\n\n",
  ps$dispersion_ranef
))

## ---------------------------------------------------------------------------
## 3. lmer reference fit
## ---------------------------------------------------------------------------
cat("\n=== lmer reference fit ===\n\n")
fit_lmer <- lme4::lmer(form_lmer, data = dat, REML = TRUE)
print(summary(fit_lmer))

## ---------------------------------------------------------------------------
## 4. Direct calls: rLMMNormal_reg_known_vcov() with BOTH engines --
##    sim_method = "DEFAULT" (the function default): exact iid draws from the
##    closed-form joint Gaussian posterior (no Markov chain, no burn-in), and
##    sim_method = "TWO_BLOCK_GIBBS": Theorem~3-calibrated Gibbs sweeps per
##    stored draw. Both engines target the exact same posterior here (every
##    Block~2 component is dNormal(), so the joint posterior is exactly
##    Gaussian) -- they differ only in *how* they draw from it, so Section 5
##    below can compare them directly against each other and against lmer.
##
## progbar/verbose match demo("Ex_12_lmerb_BigWordClub", package = "lmebayes"):
## that demo calls lmerb() without overriding progbar/verbose, and lmerb()'s
## own formals are progbar = NULL (falsy -- no bar shown) and a hardcoded
## verbose = TRUE passed to rlmerb().
## ---------------------------------------------------------------------------
fit <- rLMMNormal_reg_known_vcov(
  n            = 10000L,
  y            = design$y,
  D            = design$Z,
  group        = grp,
  W            = design$X_hyper,
  prior_list   = prior_list,
  pfamily_list = pf,
  progbar      = FALSE,
  verbose      = TRUE
)
cat(sprintf("\nsim_method_used: %s\n", fit$sim_method_used))

fit_gibbs <- rLMMNormal_reg_known_vcov(
  n            = 10000L,
  y            = design$y,
  D            = design$Z,
  group        = grp,
  W            = design$X_hyper,
  prior_list   = prior_list,
  pfamily_list = pf,
  progbar      = FALSE,
  verbose      = FALSE,
  sim_method   = "TWO_BLOCK_GIBBS"
)
cat(sprintf("sim_method_used: %s (m_convergence = %d)\n",
            fit_gibbs$sim_method_used, fit_gibbs$m_convergence))

## ---------------------------------------------------------------------------
## 5. Block 2 fixed effects: iid vs Gibbs vs ICM mean vs lmer fixef
##    (+ uncertainty)
##
## 'iid SD'/'gibbs SD' = posterior SD of gamma_k from each engine's own draws
## (the Bayesian analogue of lmer's 'Std. Error'); NOT the Monte Carlo
## SE(mean), which is much smaller and only measures how precisely n = 10000
## draws pin down the posterior mean itself. 'diff(SE)' re-expresses the iid-
## mean vs lmer gap in units of lmer's own Std. Error -- the right scale to
## judge it on (see Ex_11's Section 6 for the same table against a glmmTMB
## reference instead of lmer).
## ---------------------------------------------------------------------------
re_names <- design$re_coef_names
n_draws  <- nrow(fit$fixef[[re_names[1L]]])

fe_ref <- lmebayesCore:::.lmebayes_reference_fixef(fit_lmer)
se_ref <- sqrt(diag(lmebayesCore:::.lmebayes_reference_vcov(fit_lmer)))

## Build every row first (no printing) so the header/dashes/rows below print
## as one contiguous block -- otherwise demo()'s echo = TRUE interleaves this
## loop's own source between the header and the first data row.
rows_fe <- character(0L)
for (k in re_names) {
  dm_iid   <- colMeans(fit$fixef[[k]])
  sd_iid   <- apply(fit$fixef[[k]], 2L, sd)
  dm_gibbs <- colMeans(fit_gibbs$fixef[[k]])
  sd_gibbs <- apply(fit_gibbs$fixef[[k]], 2L, sd)
  icm_k    <- fit$fixef.mode[[k]]
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
  "\n=== Block 2 fixed effects: iid vs Gibbs vs ICM mean vs lmer fixef (+ uncertainty) ===\n\n",
  sprintf("  %-18s  %-28s  %10s  %8s  %10s  %8s  %10s  %10s  %8s  %8s\n",
          "RE component", "parameter", "iid mean", "iid SD", "gibbs mean", "gibbs SD",
          "ICM mean", "lmer", "lmer SE", "diff(SE)"),
  sprintf("  %-18s  %-28s  %10s  %8s  %10s  %8s  %10s  %10s  %8s  %8s\n",
          strrep("-", 18L), strrep("-", 28L), strrep("-", 10L), strrep("-", 8L),
          strrep("-", 10L), strrep("-", 8L), strrep("-", 10L), strrep("-", 10L),
          strrep("-", 8L), strrep("-", 8L)),
  rows_fe,
  "\n  diff(SE) = (iid mean - lmer estimate) / lmer Std. Error -- |diff(SE)| < ~1-2\n",
  "  is well within lmer's own uncertainty for that coefficient, not a discrepancy.\n",
  sep = ""
)

## ---------------------------------------------------------------------------
## 6. Random effects: MCMC mean (per group, per draw average) vs exact ICM
##
## fit$coefficients: long data.frame with one row per (draw, group), columns
## 'draw', the group-name column, and one column per RE component -- these
## are beta_j draws (the full, non-centered coefficient; see ?rLMM_reg's
## "Model and notation" section), directly comparable to lme4::coef(), not
## lme4::ranef(). fit$ranef.mode is the matching exact-ICM posterior mean.
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
n_flagged <- 0L
rows_re <- character(0L)
for (lev in grp_levs) {
  lev_chr <- as.character(lev)
  for (k in re_names) {
    mcmc_m <- re_draws_mean[[lev_chr]][[k]]
    mcmc_s <- re_draws_sd[[lev_chr]][[k]]
    icm_m  <- fit$ranef.mode[lev_chr, k]
    se_val <- mcmc_s / sqrt(n_draws)
    z_val  <- (mcmc_m - icm_m) / se_val
    flag   <- if (abs(z_val) > 3) " *" else "  "
    if (abs(z_val) > 3) n_flagged <- n_flagged + 1L
    rows_re <- c(rows_re, sprintf(
      "  %-6s  %-18s  %10.4f  %10.4f  %10.4f  %6.2f%s\n",
      lev_chr, k, mcmc_m, icm_m, se_val, z_val, flag
    ))
  }
}

total_tests <- length(grp_levs) * length(re_names)

## Single cat() call for the whole block (title + header + dashes + rows +
## footnote): demo()'s echo = TRUE prints one "> ..." per top-level R
## statement, so splitting the table across multiple cat() calls -- even
## adjacent ones -- still interleaves an echoed statement between each
## piece. One call in, one block out.
cat(
  "\n=== Random effects: MCMC mean vs exact ICM posterior mean ===\n\n",
  sprintf("  %-6s  %-18s  %10s  %10s  %10s  %6s\n",
          "group", "RE component", "MCMC mean", "ICM mean", "SE(mean)", "z"),
  sprintf("  %-6s  %-18s  %10s  %10s  %10s  %6s\n",
          strrep("-", 6L), strrep("-", 18L),
          strrep("-", 10L), strrep("-", 10L), strrep("-", 10L), strrep("-", 6L)),
  rows_re,
  sprintf(
    "\n  %d of %d tests flagged |z| > 3  (expected ~%.1f by chance at 0.3%% level)\n",
    n_flagged, total_tests, total_tests * 0.003
  ),
  "  (* |z| > 3: MCMC mean inconsistent with exact ICM posterior mean)\n",
  sep = ""
)

## ---------------------------------------------------------------------------
## 7. Random effects: lmer "full" coefficient vs sampler ranef.mode
##
## fit$ranef.mode is beta_j = W_j %*% gamma + u_j (see ?rLMM_reg's "Model and
## notation"). lme4::coef()[[k]] is NOT beta_j whenever RE component k's
## hyper-design W_k has more than an intercept column: coef() only adds
## fixef(k) + ranef(k)_j, and silently leaves any cross-level covariate's
## contribution (e.g. "(Intercept)"'s private_school/title1/free_reduced_lunch
## main effects, or "free_reduced_lunch:distracted_a1") in its own,
## group-constant column instead of folding it back into k's column. So
## coef()[[k]] alone is missing exactly the part that varies the random
## coefficient across groups by W_k's non-intercept columns -- comparing it
## directly to ranef.mode is an apples-to-oranges comparison for any RE
## component with covariates in its hyper-design (here "(Intercept)" and
## "distracted_a1"; "distracted_ppvt" is intercept-only, so unaffected).
##
## Fix (mirrors demo("Ex_24_lmerb_dGamma_BigWordClub", package = "lmebayes")):
## reconstruct the same beta_j from lmer's own fit -- ranef(k)_j (coef() minus
## its own anchor term) plus build_mu_all()'s W_j %*% gamma reconstruction,
## using lmer's fixef mapped onto W_k's columns (main effect or interaction,
## as appropriate) instead of the sampler's gamma.
## ---------------------------------------------------------------------------
cat("\n=== Random effects: lmer 'full' coefficient vs sampler ranef.mode (all groups) ===\n\n")

fe_lmer  <- lme4::fixef(fit_lmer)
coef_raw <- as.data.frame(coef(fit_lmer)[[grp_col]])
if (!identical(names(coef_raw), re_names) && all(re_names %in% names(coef_raw))) {
  coef_raw <- coef_raw[, re_names, drop = FALSE]
}

## Map each hyper-design column of W[[k]] onto the matching lmer fixef name:
## col == "(Intercept)" -> lmer's main-effect term for k itself; k ==
## "(Intercept)" -> lmer's main-effect term for col; otherwise the
## col:k / k:col interaction term.
.fe_name_for_lmer <- function(k, col, fe) {
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

fixef_lmer <- lapply(re_names, function(k) {
  cols_k <- colnames(design$X_hyper[[k]])
  fe_nms <- vapply(cols_k, .fe_name_for_lmer, character(1L), k = k, fe = fe_lmer)
  miss <- is.na(fe_nms) | !fe_nms %in% names(fe_lmer)
  if (any(miss)) {
    stop(
      "lmer fixef missing term(s) for W[[", k, "]]: ",
      paste(cols_k[miss], collapse = ", "),
      call. = FALSE
    )
  }
  mu_k <- vapply(fe_nms, function(nm) unname(fe_lmer[nm]), numeric(1L))
  names(mu_k) <- cols_k
  mu_k
})
names(fixef_lmer) <- re_names

## coef_anchor[k]: the fixef term coef()[[k]] already includes (fixef(k)
## itself), so it must be subtracted before adding the full W_j %*% gamma_lmer
## reconstruction below (otherwise k's own intercept contribution is double
## counted).
coef_anchor <- vapply(re_names, function(k) {
  if (k == "(Intercept)") unname(fe_lmer["(Intercept)"]) else unname(fe_lmer[k])
}, numeric(1L))

mu_all_lmer <- build_mu_all(design, fixef_lmer)$mu_all

cmp <- do.call(rbind, lapply(grp_levs, function(lev) {
  data.frame(
    group       = lev,
    re_coef     = re_names,
    lmer_full   = vapply(re_names, function(k) {
      mu_all_lmer[k, lev] + (coef_raw[lev, k] - coef_anchor[[k]])
    }, numeric(1L)),
    sampler_icm = unname(fit$ranef.mode[lev, re_names]),
    stringsAsFactors = FALSE
  )
}))
cat(
  "  lmer_full = build_mu_all(design, fixef_lmer)$mu_all + ",
  "(coef(fit_lmer) - coef_anchor); compare to sampler_icm (fit$ranef.mode):\n\n"
)
print(round(cmp[, c("lmer_full", "sampler_icm")], 4), row.names = paste(cmp$group, cmp$re_coef, sep = "::"))

## ---------------------------------------------------------------------------
## 8. Cross-check: sim_method = "DEFAULT" (exact iid) vs "TWO_BLOCK_GIBBS"
##    random effects (fixed effects are already cross-checked in Section 5)
##
## fit/fit_gibbs's ranef.mode should agree closely -- this is exactly the
## check used to catch and verify the fix for the ranef.mode bug described
## above (see NEWS.md): before the fix, the iid route's ranef.mode was the
## *last draw* (noisy), not a point estimate, and disagreed with
## TWO_BLOCK_GIBBS's ICM-mode ranef.mode by tens of units; after the fix,
## both engines' ranef.mode match to within Monte Carlo error.
## ---------------------------------------------------------------------------
cat("\n=== Cross-check: sim_method = \"DEFAULT\" vs \"TWO_BLOCK_GIBBS\" (random effects) ===\n\n")
cat("-- Random effects: ranef.mode agreement (all groups) --\n\n")
ranef_cmp <- do.call(rbind, lapply(grp_levs, function(lev) {
  data.frame(
    group        = lev,
    re_coef      = re_names,
    iid_mode     = unname(fit$ranef.mode[lev, re_names]),
    gibbs_mode   = unname(fit_gibbs$ranef.mode[lev, re_names]),
    stringsAsFactors = FALSE
  )
}))
print(
  round(ranef_cmp[, c("iid_mode", "gibbs_mode")], 4),
  row.names = paste(ranef_cmp$group, ranef_cmp$re_coef, sep = "::")
)
