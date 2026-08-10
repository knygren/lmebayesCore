# lmebayesCore (development version)

* **`rlmerb()` / `rglmerb()` prior formals align with `rLMM*` /
  `rGLMM*`.** Replaced packed `prior` + `group.dispersion` with
  `pfamily_list` and `dispprior_list` (same shapes as the matrix engines;
  bare scalars / `dGamma()` / `dGamma_list()` still accepted). Also added
  `offset = NULL` / `weights = 1` (same defaults as `rLMM*` / `rGLMM*`;
  omitted args still inherit from `design` when present). Formula drivers
  in **lmebayes** updated accordingly.

* **Further Core helpers no longer exported** (help kept as
  `\keyword{internal}`; call via `:::`): `normalize_group()`,
  `build_mu_all()`, `matrix_args_lmm()`, `lmerb_posterior_mean()`,
  `glmerb_posterior_mode()`, `lmerb_posterior_covariance()`,
  `rGLMM_sweep()`, `rGLMM_Re_Draw()`, plus the earlier two-block rate /
  pilot / Block~2 helpers. Downstream **lmebayes** paths updated to
  `:::` where needed.

* **Two-block helpers no longer exported.** Documented engines and
  calibration helpers (`two_block_rNormal_reg`, `two_block_rate*`,
  `two_block_tv_bound` / `two_block_l_for_tv`, pilot-cost stack,
  `two_block_align_b_to_xhyper*`, `two_block_block2_one_chain*`, and
  dotted aliases) stay in help as internals; call via `:::` if needed.
  Public entry remains `rlmerb` / `rglmerb` / `rLMM*` / `rGLMM*`.
  S3 `print` methods for rate objects stay registered.

* **Observation-partition APIs use `group`, not `block`.** Renamed
  `rNormal_reg_block` / `rNormalGLM_reg_block` →
  `rNormal_reg_group` / `rNormalGLM_reg_group` (formal `group=`;
  return fields `group_info` / `group_results`). Also
  `normalize_block` → `normalize_group`, `Prior_SetupBlock` →
  `Prior_SetupGroup`, doc topic `simfuncs_group`. Gibbs engine names
  `two_block_*` stay (parameter blocks). Compiled
  `.block_rNormal*_cpp` entry points keep their old names.

* **`Prior_Setup_GLMM` marginal dispersion centers.** Pooled Block~1
  \(\sigma^2\) (`dispformula = ~1`) and Block~2 \(\tau^2_k\) Gamma priors
  are now centered on A12 §3.3.4 marginal RSS via
  `compute_gaussian_prior` (same idea as per-group measurement /
  Gaussian `Prior_Setup`), not raw lmer residual / VarCorr plug-ins.
  Pooled \(\hat\sigma^2\) aggregates per-group \(S_{\mathrm{marg},j}\)
  (Part VI \(\Omega_j\) included); \(\tau^2_k\) uses the hyper-regression
  \(b_{\cdot k}\sim N(W_k\gamma_k,\tau^2_k I)\). Classical plugs remain in
  `attr(group.dispersion,"classical")` and `pop.dispersion.ref`.

* **`Prior_Setup_GLMM`: drop `group.dispersion=` pin.** Removed the
  argument that overrode / mean-matched Gamma centers to a user plug-in
  (pooled or per-group). Centers always come from the marginal
  calibration above; pass `ps$group.dispersion` into
  `rlmerb()`/`rglmerb()` only when you want a fixed-\(\sigma^2\) sampling
  route.

* **Configure policy:** Removed `tools/rcpp_include.R` /
  `tools/patch_rcpp_function_h.R` and the
  `glmbayes_getRegisteredNamespace` compatibility shim
  (`src/glmbayes_getRegisteredNamespace.{cpp,h}`) from `configure` and
  `configure.win`. Builds now rely on standard `LinkingTo: Rcpp` and
  `Rcpp (>= 1.1.1)` instead of probing/patching `Rcpp/Function.h`.
  Synced `.github/workflows/rhub.yaml` with **glmbayesCore** (drops dead
  patch/GitHub-Rcpp matrix logic).

* **`rLMMindepNormalGamma_reg_*_v2` stubs + rate derivation.** New exports
  `rLMMindepNormalGamma_reg_known_vcov_v2()` and
  `rLMMindepNormalGamma_reg_estimated_vcov_v2()` (same formals as v1) for a
  future Gibbs partition that updates observation \(\sigma^2\) with the
  population block instead of jointly with group coefficients. Currently
  `stop()` with a pointer to
  `inst/DERIVATION_sigma2_with_block2_v2.md` (lambda*/eigenvalue changes
  relative to §10/§14 of `BLOCK_GIBBS_ERGODICITY_ING.md`). Helpers
  `.lmebayes_ing_measurement_prior_list_v2()` /
  `_group_v2()` drop unused ING `mu`/`Sigma`.

* **`rLMM*` / `rGLMM*`: `offset` / `weights` formals (Phase 1 echo).**
  Matrix engines now accept glmbayes-style `offset = NULL` and
  `weights = 1`, normalize them to `length(y)`, and return
  `$prior.weights`, `$offset`, and `$offset2` (also stored on
  `$design`). `rlmerb` / `rglmerb` forward `design$weights` /
  `design$offset` from `model_setup`. **Not yet used** in ICM or
  Gibbs sweeps (sampling still assumes unit weights and zero offset).

* **`rLMM*` / `rGLMM*` return objects use group/pop names.** Public slots
  rename from lme4-style `coefficients`/`fixef`/`ranef.*` to package
  notation: `groupef` / `groupef.mode` / `groupef.iters` (group
  \(\beta_j\)) and `popef` / `popef.mode` / `popef.init` /
  `popef.dispersion` / `popef.iters` (population \(\gamma\)). Observation
  residual variance is `$group.dispersion` / `$group.dispersion.mean`
  (replacing `$dispersion_ranef` and `$sigma2`). Pilot metadata nests
  under `$pilot`; engine labels and ICM info nest under
  `$convergence_info`. Unused slots dropped. Shared assembler
  `.lmebayes_assemble_reg_result()` enforces order.

* **`rLMM*` / `rGLMM*` engines: `prior_list` renamed to `dispprior_list`**
  and moved after `pfamily_list` in the formals
  (`…, W, pfamily_list, dispprior_list, …`). Return objects store the
  Block~1 measurement prior as `$dispprior_list`. Call sites
  (`matrix_args_lmm()`, `.lmebayes_matrix_args_glmm()`, demos, tests)
  updated. `\value` documentation added for `?rLMM_reg` and
  `?rGLMM_reg`.

* **`pfamily_list()` generic now includes `ptypes`**
  (`pfamily_list(object, ptypes = "dNormal", ...)`), so the prior-family
  type(s) are part of the shared interface rather than a method-only
  formal passed through `...`.

* **`pfamily_list()` returns class `"pfamily_list"`** with
  `print.pfamily_list()` that prints each component via
  `print.pfamily`, one after another. Use
  `print(x, components = ...)` to print a subset of population /
  group-effect models (e.g. `"(Intercept)"`, slopes); `NULL` prints
  all.

* **`dGamma_list()` returns class `"dGamma_list"`** with
  `print.dGamma_list()` that prints each group's `dGamma()` via
  `print.pfamily` only (does not dump list attributes such as
  `window_diagnostics` or the stored `glmmTMB` fit). Use
  `print(x, groups = ...)` to print a subset of group names (or
  indices); `NULL` prints all.

* **`print.model_setup()` section headers relabeled** to Level 1
  (Within-Group) / Level 2 (Across-Group) / Level 2 Rank, with
  `Group-effect coefficients` instead of `RE predictors` /
  Measurement Model / Random Effects Model wording.

* **`print.Prior_Setup_GLMM()` no longer dumps the per-group Block~1
  `sigma^2` calibration table** (or the `group.pwt_calibration` table).
  Those details stay on `group.ing_prior` / `dGamma_list()`; the print
  keeps a one-line group count and pwt / `max_disp_perc` summary.
  Also removed the construction-time A12 rate_gamma vs rate console dump.
  Print layout now uses spaced sections (Setup / Group dispersion /
  Design check / group.Sigma / Population priors) with a short description
  under each header (no Block 1/2 labels in the printout). Group-dispersion
  prior detail is omitted from the main print; it points to
  `group.ing_prior` / `print(dGamma_list(...))` instead.

* **`glmmTMB` moved from Suggests to Imports** (used from `R/` for
  per-group dispersion reference fits). **`bayesrules` stays in
  Suggests**: it is only needed for example/demo data
  (`big_word_club`). Help examples that use it are wrapped in
  `requireNamespace("bayesrules", quietly = TRUE)` so
  `example()` skips quietly when the Suggests package is not installed.

* **Rename: `Prior_Setup_lmebayes()` is now `Prior_Setup_GLMM()`.**
  The returned S3 class is `"Prior_Setup_GLMM"` (was
  `"lmebayes_prior_setup"`). S3 methods
  `pfamily_list.Prior_Setup_GLMM()`, `dGamma_list.Prior_Setup_GLMM()`,
  and `print.Prior_Setup_GLMM()` follow the new class name. Call sites
  and re-exports (including in \pkg{lmebayes}) must use the new name.

* **New: `two_block_rate_ing()`, an extended (lambda, Omega)-aware LOCAL
  rate diagnostic for the two-block Gibbs chain.**
    - Extends `two_block_rate()`'s `(gamma, beta)`-only Hessian with the
      diagonal RE-precision (`lambda_p = 1/tau_p^2`) and/or measurement-
      precision (`Omega`/`Omega_j = 1/sigma^2`/`1/sigma_j^2`) blocks derived
      in `inst/BLOCK_GIBBS_ERGODICITY_ING.md`, evaluated at a single
      reference state (residuals `u`/`e` and precision values supplied via
      the new `lambda_ing`/`omega_ing` arguments).
    - Diagnostic only: this is a state-dependent local rate, not a
      certified total-variation bound (see `?two_block_rate_ing` and
      `print.two_block_rate_ing()`'s printed caveat) -- it does not change
      any sampler's calibrated `m_convergence` and is not called by any
      exported sampler.
    - New internal helpers `.lmebayes_reference_u()` and
      `.lmebayes_reference_group_residuals()` extract the RE-level
      (`u_jp`) and data-level (`e_j`) residual plug-ins from an `lmer`/
      `glmmTMB` reference fit.
    - Wired into `demo("Ex_12_...")` (lambda-extension), `demo("Ex_13_...")`
      (Omega-extension), and `demo("Ex_14_...")` (combined) as a reported
      diagnostic printed next to the existing certified corner-bound rate.
    - Fix: unlike `two_block_rate()`'s exact `(gamma, beta)` system (where
      Remark~8 guarantees `lambda_star < 1` everywhere, so hitting the
      numerical ceiling means a computation bug), the extended
      `(beta, lambda/Omega)` Hessian is only guaranteed positive definite
      *near its own joint mode* -- evaluating it at a reference state where
      a group's residual is inconsistent with its plugged-in precision (e.g.
      an external reference fit that poorly predicts that group's data) can
      legitimately produce `lambda_star >= 1`. `two_block_rate_ing()` now
      reports this (with a `warning()`) instead of erroring
      (`.two_block_gen_eigen()` gained a `strict` argument).
    - New internal helpers `.lmebayes_posterior_mean_group_coef()`,
      `.lmebayes_posterior_group_residuals()`, and `.lmebayes_posterior_u()`
      evaluate `lambda_ing`/`omega_ing` at a *completed sampler's own*
      posterior-mean state (`fit$coefficients`/`fit$fixef`/
      `fit$dispersion_ranef.mean`/`fit$fixef.dispersion`) instead of an
      external `lmer`/`glmmTMB` reference fit, so the plugged-in precision
      and its paired residual always come from the same mutually-consistent
      state; `demo("Ex_13_...")` and `demo("Ex_14_...")` now use these.
    - New `demo("Ex_13b_...")`: a copy of `demo("Ex_13_...")` using the
      (unimplemented-in-`dGamma_list()`) Part~VI extension of
      `inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md` -- integrating out the
      per-group prior mean `mu_j` via a model-derived
      `Omega_j = W_j Sigma_fixef W_j'` -- for comparison against `Ex_13`'s
      production calibration.
    - New internal helpers `.lmebayes_fit_at_draw()` and
      `.two_block_rate_ing_over_draws()` evaluate the extended
      `two_block_rate_ing()` diagnostic at *every* main-stage draw already
      returned by a completed fit (rather than a single reference state),
      treating the draws as approximate posterior samples and tracking the
      empirical worst-case `lambda_star`/eigenvalues across them -- the
      main-stage analogue of `.two_block_pilot_ub_from_coefficients()`'s
      pilot-draw `pmax()` scan, applied to the extended `(lambda, Omega)`
      system. `demo("Ex_12_...")`, `demo("Ex_13_...")`, `demo("Ex_13b_...")`,
      and `demo("Ex_14_...")` now print this alongside the existing
      single-reference-state diagnostic. Diagnostic only; does not feed
      back into any sampler's calibration.

* **New: `model_setup()` gains a `dispformula` argument and additive
  `design$glmmTMB_fit` field.**
    - `dispformula = ~1` (default, pooled) preserves current behavior
      exactly. `dispformula = ~<group_name>` (matching the random-effects
      grouping factor) additionally fits a `glmmTMB::glmmTMB()` per-group-
      dispersion reference model (Gaussian models only; errors otherwise)
      and stores it as a **new, additive** `design$glmmTMB_fit` field.
      `design$lmer_fit`/`design$glmer_fit` are never touched by this --
      they stay the plain pooled-dispersion fit in every case, exactly as
      documented before.
    - `Prior_Setup_GLMM()` now passes its own `dispformula` argument
      through to this new `model_setup()` argument and reuses
      `design$glmmTMB_fit` as its `fit_ref`/`dispersion_fit` calibration
      reference, instead of independently fitting `glmmTMB` a second time.
    - Internal: `.Prior_Setup_GLMM_dispformula_kind()` moved from
      `Prior_Setup_GLMM.R` to `glmmtmb_reference_helpers.R` and renamed
      to `.lmebayes_dispformula_kind()` so `model_setup()` can share the
      same `dispformula` classification/validation logic; behavior is
      unchanged. New internal dispatcher
      `.lmebayes_extract_reference_variance_components()` picks
      `extract_mer_variance_components()` or
      `extract_glmmtmb_variance_components()` based on the fit's class.

* **Breaking: `model_setup()`'s (and `extract_re_hyper_matrices()`'s)
  returned `"model_setup"` object renames its design-matrix fields to match
  the `rLMM_reg()`/`rGLMM_reg()`/`check_identifiability()` `D`/`W`/`group`
  convention.**
    - `Z` -> `D` (level-1 random-effect design matrix), `X_hyper` -> `W`
      (named list of level-2 hyper-design matrices), `groups` -> `group`
      (grouping factor). `re_coef_names` and `group_name` are unchanged.
    - This is a clean rename with no deprecated alias: code that read
      `design$Z`, `design$X_hyper`, or `design$groups` (from `model_setup()`,
      `extract_re_hyper_matrices()`, or the matching internal "design" lists
      built by `rLMM_reg()`/`rGLMM_reg()`/`rLMMNormal_joint_iid()`/staged
      sweep engines) must switch to `design$D`, `design$W`, `design$group`.
    - Downstream consumers updated to match: `matrix_args_lmm()`/
      `matrix_args_glmm()`, `priors_from_pfamily_list()`, `build_mu_all()`
      (and its `.lmerb_validate_design()` field check), `lmerb_posterior_mean()`
      / `glmerb_posterior_mode()` / `lmerb_posterior_covariance()`,
      `two_block_block2_one_chain()` / `two_block_block2_one_chain_cpp()`,
      `Prior_Setup_GLMM()`, `dGamma_list()`, `pfamily_list()`,
      `rlmerb()`/`rglmerb()`, and the `plot_var_convergence()`/
      `plot_mean_convergence()` fit-object methods. The C++ boundary in
      `src/two_block_block1.cpp` (`design["Z"]`/`design["groups"]`/
      `design["X_hyper"]`) was updated to match, though its only current
      caller in R (`.two_block_block1_one_chain_cpp()`) is not on any live
      call path.

* **New: `check_identifiability()`, exported standalone Level-1/Level-2
  identifiability and estimability check.**
    - Extracted from `model_setup()`'s inline rank/estimability block (the
      "two-step identifiability assessment": per-group `D_j` rank and
      estimability, then per-RE-coefficient `W[[k]]` rank restricted to
      estimable groups). `model_setup()` now calls `check_identifiability()`
      and copies its fields onto the returned design object (`D = design$D`,
      `group = design$group`, `W = design$W`); its behavior is unchanged.
    - Argument names (`D`, `group`, `W`) match the recently-renamed
      `rLMM_reg()`/`rGLMM_reg()` matrix-level conventions -- `y`/`D`/`group`/
      `W` share their `@param` documentation with `rLMM_reg()` via
      `@inheritParams`, and `W[[k]]` uses the same row-`j`-\eqn{\leftrightarrow}
      -`levels(group)[j]` positional convention (no `rownames()` requirement).
    - Runs directly on hand-built `(y, D, group, W)` without a `lmer`/
      `glmer` reference fit, so matrix-level `rLMM_reg()`/`rGLMM_reg()`
      callers can check identifiability before sampling.
    - Validates its inputs strictly: `D` must have unique, non-empty
      `colnames(D)` (there is no separate `re_coef_names` argument); `group`
      must be a `factor` (no separate `group_levels` argument); `W` must be
      a named list whose names match `colnames(D)` exactly (as sets) and
      whose elements each have `length(levels(group))` rows; `group_name`
      must be resolvable from the argument or `attr(group, "group_name")`.
      Mismatches raise an informative error immediately.
    - Returns `re_rank`, `re_estimable`, `re_glm_check`, `hyper_rank`,
      `hyper_deficient`, `rank_ok`, plus `group_levels` and `group_name`.
    - `.lmebayes_block_glm_estimable()` and
      `.lmebayes_glm_estimable_precheck()` moved from `model_setup.R` to the
      new `R/check_identifiability.R` (internal; no behavior change).

* **`model_setup()` now performs the binomial classical-glm MLE-existence
  check itself, and the Level~2 hyper-design rank check is restricted to
  *estimable* groups (not merely full-rank ones).**
    - `model_setup()`'s return value gains `re_estimable` (named logical,
      per group) and `re_glm_check` (per-group diagnostic data frame, or
      `NULL` for non-binomial families). For `family = binomial()`, a group
      that is algebraically full-rank (`re_rank`) but has no finite
      classical-`glm(y ~ Z_j - 1, family = binomial)` MLE (complete or
      quasi-complete separation) is now flagged `re_estimable = FALSE`; other
      families currently set `re_estimable` equal to `re_rank` (no glm check
      yet -- Poisson and other families are a planned follow-up).
    - This check previously lived inside `Prior_Setup_GLMM()`
      (`.lmebayes_block_glm_estimable()`, attached post-hoc to the `design`
      object `model_setup()` had already returned). It has moved into
      `model_setup()` itself; `Prior_Setup_GLMM()` now just reads
      `design$re_estimable`/`design$re_glm_check` instead of recomputing
      them. No argument was added to opt out -- for `family = binomial()`
      this check always runs.
    - **Behavior change:** `model_setup()`'s Level~2 `hyper_rank` check (is
      `X_hyper[[k]]` full column rank?) now restricts to `re_estimable`
      groups instead of `re_rank` groups. For `family = binomial()` with any
      full-rank-but-non-estimable (separated) group, this can change
      `hyper_rank`/`hyper_deficient`/`rank_ok` relative to previous versions,
      since such groups no longer count toward Level~2 identifiability (a
      group without a finite classical MLE supplies no real information
      about its random-effect coefficients, so it shouldn't count as
      identifying the level-2 hyperparameters either). For non-binomial
      families `re_estimable == re_rank`, so `hyper_rank` is unchanged.
    - `print.model_setup()` now reports the glm-MLE count and any
      full-rank-but-not-estimable groups for `family = binomial()`, and its
      "Hyper-Design Rank" section now labels the restriction set "estimable"
      rather than "full-rank".

* **The per-group estimability check (`re_estimable`/`re_glm_check`) now
  also covers `poisson()` and `Gamma()`, and `gaussian()` gains a
  residual-degrees-of-freedom check.**
    - `poisson()`: an algebraically full-rank group with an all-zero count
      response is flagged `re_estimable = FALSE` (log-link intercept MLE
      diverges to `-Inf`), otherwise the usual classical-`glm`
      fit/coefficient/`vcov` finiteness checks apply (as for `binomial()`).
    - `Gamma()`: a group with any non-positive response value is flagged
      `re_estimable = FALSE` (invalid domain for `Gamma()`); a group with
      `n_j <= p_j` (no residual degrees of freedom) is also flagged
      `re_estimable = FALSE`, since the dispersion/shape parameter would not
      be estimable even though the mean-model MLE exists.
    - `gaussian()`: **behavior change.** Previously `re_estimable` was set
      equal to `re_rank` unconditionally. A classical OLS coefficient MLE
      always exists given full column rank, but the group-level residual
      *dispersion* additionally requires `n_j > p_j`; a full-rank group with
      `n_j <= p_j` (a perfect/saturated fit) is now flagged
      `re_estimable = FALSE`. No `lm()`/`glm()` fit is attempted for
      `gaussian()` -- this is a cheap arithmetic check on already-known
      `n_j`/`p_j`. `re_glm_check` is therefore no longer always `NULL` for
      `gaussian()`.
    - This can change `hyper_rank`/`hyper_deficient`/`rank_ok` for
      `Gamma()`/`gaussian()` models with small per-group sample sizes
      relative to the number of random-effect predictors, since
      non-estimable groups no longer count toward Level~2 identifiability.
    - `print.model_setup()`/`print.Prior_Setup_GLMM()`'s "Full-rank with
      glm MLE" line is renamed "Full-rank & estimable" and now fires for any
      family with a non-`NULL` `re_glm_check` (not just `binomial()`).
    - `quasipoisson()`/`quasibinomial()` and other families are unchanged
      (`re_estimable == re_rank`, `re_glm_check = NULL`).

* **New: `plot_mean_convergence()`, and fit-object (S3) methods for
  `plot_var_convergence()`/`plot_mean_convergence()`/`plot_sweep_history_diag()`.**
    - `plot_mean_convergence()` is the mean-bias companion to
      `plot_var_convergence()`, implementing Claim~1 of the two-block Gibbs
      ergodicity reference (`inst/BLOCK_GIBBS_ERGODICITY.md`): after `l`
      inner sweeps, each chain's mean is biased toward its start by
      `A^l`, so `(mean(l) - mean_ref) / sd(l)` shrinks toward 0 as `l`
      grows. Same two modes as `plot_var_convergence()`
      (named-coefficient `coef_focus` traces, or `whitened = TRUE`
      eigenvalue traces of the whitened mean-deviation vector), the same
      automatic exact-vs-empirical reference resolution (`lmerb_posterior_mean()`
      when `design`/`measurement_prior_list` are supplied and both dispersion
      and the RE vcov are fixed; the empirical last-sweep cross-chain mean
      otherwise), and an analogous `n_chains`/`conf_level` confidence band
      (Gaussian/`qt()`-based, vs. `plot_var_convergence()`'s chi-squared/`qf()`).
    - `plot_var_convergence()`, `plot_mean_convergence()`, and
      `plot_sweep_history_diag()` are now S3 generics. The previous
      `hist`-first-argument behavior is preserved as each generic's
      `.default` method (unchanged signature/behavior). New fit-object
      methods let you call any of the three directly on a fitted object --
      `rlmerb()`, `rglmerb()`, `rLMMNormal_reg_known_vcov()`/
      `rLMMNormal_reg_estimated_vcov()`, `rLMMindepNormalGamma_reg_known_vcov()`/
      `rLMMindepNormalGamma_reg_estimated_vcov()`, and
      `rGLMM_reg_known_vcov()`/`rGLMM_reg_estimated_vcov()` -- instead of
      manually extracting `fit$sweep_history`/`fit$design` and building a
      `measurement_prior_list` by hand. A new `stage = c("main", "pilot")`
      argument selects `fit$sweep_history` or `fit$pilot$sweep_history`;
      `n_chains` defaults to the *correct* count for that stage (`fit$n` for
      `"main"`, `fit$pilot_chisq$n_pilot` for `"pilot"` -- these can differ,
      e.g. when `n_pilot` is calibrated from `gap_tol`). The exact reference
      (vs. the empirical `Var_final`/`mean_final` fallback) is resolved
      automatically per fit: available only for Gaussian fits with *both*
      dispersion and the Block~2 vcov fixed (`rLMMNormal_reg_known_vcov()`
      and `rlmerb()`/`rglmerb()` routed there); never available for
      `rLMMindepNormalGamma_reg_*()` (dispersion always estimated per group)
      or `rGLMM_reg_*()` (non-Gaussian in every real `glmerb()` call path).
    - Bug fix uncovered while wiring this up: `rLMMNormal_reg_known_vcov()`'s
      exact-iid engine (`sim_method = "DEFAULT"`, the function's actual
      default) was missing `$design` on the returned fit -- only the
      `sim_method = "TWO_BLOCK_GIBBS"` route had it. Fixed (`.rLMMNormal_reg_run_iid()`
      now attaches the same `design` list the Gibbs route does), though this
      route still has no `sweep_history` (no Gibbs sweeps to plot at all;
      the new fit-object plotting methods now error with a clear message
      pointing at `sim_method = "TWO_BLOCK_GIBBS"` instead of failing on a
      missing field).
    - Demo simplification: `Ex_10`/`Ex_11` (`rLMMNormal_reg_known_vcov()`,
      the fully-known case) now call `plot_mean_convergence(fit_gibbs, ...)`/
      `plot_var_convergence(fit_gibbs, ...)` directly instead of building
      `measurement_prior_list` by hand via `.two_block_measurement_prior_list()`/
      `.rLMM_P_from_pfamily_list()`. `Ex_12`/`Ex_13`/`Ex_14` (estimated
      vcov and/or estimated dispersion) now loop
      `for (stg in c("pilot", "main")) plot_mean_convergence(fit, stage = stg); plot_var_convergence(fit, stage = stg)`
      instead of a hand-built `list(fit$pilot$sweep_history, fit$sweep_history)`
      -- as a side effect this also fixes `Ex_12`'s pilot-stage confidence
      band, which had hardcoded the main-stage `n_chains = 3000` for both
      stages even though its pilot stage actually runs a different,
      `gap_tol`-calibrated chain count.
    - New `split`/`max_whitened` arguments on `plot_var_convergence()`/
      `plot_mean_convergence()` (and forwarded through `...` on every
      fit-object method above). `split = "auto"` (the new default) draws
      one *separate* chart per group -- its own plot page/window, never
      combined via `layout()`/`par(mfrow = ...)` into one page -- instead
      of a single chart with every coefficient/eigenvalue overlaid:
      for named-coefficient mode, one chart for the intercept's own
      hyper-predictors ("Intercept predictors") and one for every other
      random-effect component's hyper-predictors combined ("Slope
      predictors"); for `whitened = TRUE`, consecutive batches of at most
      `max_whitened` eigenvalues each (default `4`). `split = "none"`
      restores the previous single-combined-chart behavior.

* **New: `component = "precision"` on `plot_var_convergence()`/
  `plot_mean_convergence()`, plotting the estimated random-effects
  *precision* (`1/tau^2_k`, the diagonal of `Sigma_ranef^-1`) instead of
  the Block~2 fixed-effects/hyperparameter series.**
    - `rGLMM_sweep()` and `.rGLMM_sweep_ing_block1()` (the estimated-vcov
      engines behind `rLMMNormal_reg_estimated_vcov()`,
      `rLMMindepNormalGamma_reg_estimated_vcov()`, and
      `rGLMM_reg_estimated_vcov()`) now also snapshot the cross-chain
      mean/SD of `1/tau2_k` after every inner sweep, for every ING
      (`dIndependent_Normal_Gamma`) Block~2 component (`dNormal` components
      have a fixed dispersion, so there is nothing to track for them).
      Stored as `sweep_history$disp_table`, same 5-column shape as the
      existing `table` (`re_component`, `covariate = "precision"`, `sweep`,
      `mean`, `sd`); a 0-row frame for known-vcov engines or an
      estimated-vcov fit with zero ING components. `print()` on a
      `two_block_sweep_history` now also shows a short "Block 2 RE
      precision (1/tau^2)" table when `disp_table` has rows.
    - Precision (`1/tau^2_k`), not the variance `tau^2_k` itself, is
      plotted: `E[1/tau^2_k]` stays well-defined under a weak/vague prior
      where `E[tau^2_k]` need not be (same reasoning already used by
      `.two_block_tau2_start_from_dispersion_draws()`'s `1/mean(1/x)`
      plug-in). This is unrelated to the observation-level dispersion
      estimated by `rLMMindepNormalGamma_reg_*()` (despite both being
      called "dispersion" internally) -- `component = "precision"` is
      always about the RE variance-covariance.
    - No exact reference exists yet for an estimated `Sigma_ranef`, so
      `component = "precision"` always uses the empirical last-sweep value
      as its denominator/reference (`design`/`measurement_prior_list` are
      ignored). Only the cross-chain *variance* of each component's
      precision is captured (no cross-chain covariance, either between
      precision components or between precision and fixed effects), so
      `whitened = TRUE` is not yet supported for `component = "precision"`
      and errors with a clear message; likewise `component = "precision"`
      errors clearly when `disp_table` is empty (e.g. a known-vcov fit).
    - Deferred (documented as future work in
      `inst/BLOCK_GIBBS_ERGODICITY.md`, not implemented): capturing the
      precision-with-precision and precision-with-fixed-effects cross-chain
      covariances (needed for `whitened = TRUE` support here), a joint
      Hessian/precision spanning fixed effects, random effects, *and* the
      diagonal RE precision together, and re-deriving/adjusting the
      Claim~1/Claim~3 eigenvalue bounds for a 3-block (or more) Gibbs cycle
      where the vcov itself is also sampled.

* **Bug fix: `plot_var_convergence(..., engine = "base")` could error
  with `Error in graphics::par(old_par) : invalid value specified for
  graphical parameter "pin"`** when restoring its plotting parameters on
  exit. The legend-panel `layout()` helper (`.convergence_plot_base()`)
  captured a full `graphics::par(no.readonly = TRUE)` snapshot -- including
  device-geometry-derived entries like `pin` -- and force-restored all of it
  afterward; if the plotting device's actual geometry at restore time no
  longer matched what was captured (e.g. a resized RStudio Plots pane, or
  simply a small device), restoring the old absolute-inches `pin` value
  could conflict with the device's current size and error. Only `mar` (plus
  the layout panel structure, already reset via `layout(1L)`) is ever
  changed by this function, so only `mar` is now saved/restored -- the same
  pattern `plot_sweep_history_diag()` already used for its own `par()` calls
  in this file.

* **Demo: `Ex_10_rLMM_known_dispersion_known_vcov_BigWordClub` and
  `Ex_11_rLMM_known_dispersion_vector_known_vcov_BigWordClub` now include
  `plot_mean_convergence()`/`plot_var_convergence()` (non-whitened and
  whitened) for the `sim_method = "TWO_BLOCK_GIBBS"` fit**, right after each
  demo's Block~2 fixed-effects comparison table -- enabled by the
  `rGLMM_sweep()` engine swap below, which is what makes
  `fit_gibbs$sweep_history` non-`NULL` for this route in the first place.
  Both demos are the fully-known (dispersion *and* Sigma_ranef fixed) case,
  so the exact reference mean/covariance (via `lmerb_posterior_mean()`/
  `lmerb_posterior_covariance()`) is resolved automatically by the
  fit-object plotting methods (see above).

* **Internal: `rLMMNormal_reg_known_vcov_two_bg()` (the
  `sim_method = "TWO_BLOCK_GIBBS"` route of `rLMMNormal_reg_known_vcov()`,
  including its known per-group dispersion *vector* case) now runs its
  two-block Gibbs sweeps sweeps-outer/chains-inner via `rGLMM_sweep()`**,
  the same batch engine already used by `rLMMNormal_reg_estimated_vcov()`
  and the ING routes, instead of the old chains-outer/sweeps-inner C++
  driver (`two_block_rNormal_reg()` / `two_block_rNormal_reg_v2_cpp_export()`).
  User-visible consequences:
    - `$sweep_history` is now populated on the returned object (previously
      always `NULL` for this route), enabling
      `plot_sweep_history_diag()`/`plot_var_convergence()` and
      `print()` on the sweep-by-sweep Block~2 fixef table -- mirroring the
      estimated-vcov/ING routes.
    - Progress bars are now the same nested
      `"[main] sweep m/M RE:"` / `"[main] sweep m/M fixef:"` bars used by
      the other sweep-outer engines (previously a single bar advancing
      over chains, with no per-sweep breakdown). As with those engines,
      `verbose = TRUE` alone now also enables the progress bar (previously
      `verbose` only affected ICM/calibration `cat()` messages for this
      route).
    - `fixef.mu` (`mu_all_last` internally) is now the fitted mean
      response evaluated at the **across-chains mean** of the final-sweep
      fixed effects, rather than an arbitrary single chain's (the last
      one's) final draw.
  Draws themselves are statistically unchanged (same exact-conditional
  two-block Gibbs kernel, same ICM start and Theorem~3 `m_convergence`
  calibration); this is purely an implementation/engine change, covered by
  `test-sim-method-iid.R`'s existing iid-vs-Gibbs Monte Carlo agreement
  checks. The legacy `rLMMindepNormalGamma_reg()` outer-loop driver, which
  also calls this shared pipeline internally (`n = 1` chain per outer
  draw), is unaffected in behavior.

* **Enhancement: `plot_var_convergence()` gains an optional
  `n_chains` argument (with `conf_level`, default `0.95`) that draws a pair
  of horizontal dotted reference lines**: a naive confidence band around 1
  for `Ratio(l)`, under the null that the true ratio at that sweep already
  equals 1. The band's distribution depends on whether `Var_ref` is exact
  (`design`/`measurement_prior_list` supplied) or empirical (`Var_final`):
    - Exact `Var_ref` (a fixed/known number): only the numerator `Var(l)`
      has sampling error, so
      \eqn{\mathrm{Ratio}(l) \sim \chi^2_{n_{chains}-1} / (n_{chains}-1)}.
    - Empirical `Var_ref` (itself a sample variance across the same
      `n_chains` chains, at the final sweep): both numerator and
      denominator carry sampling error, so the classic two-sample
      variance-ratio result applies instead --
      \eqn{\mathrm{Ratio}(l) \sim F(n_{chains}-1,\, n_{chains}-1)}, which is
      *wider* than the chi-squared band. Both cases are computed via a
      single `stats::qf()` call (`df2 = Inf` for the exact case, which
      reduces exactly to `qchisq(p, df1) / df1`).

  This is exact for the non-whitened (named-coefficient) series and only an
  approximation for `whitened = TRUE` (whitened eigenvalues are Wishart-,
  not marginally chi-squared/F, distributed, and tend to spread out more
  than this band implies); also a further approximation in the empirical
  case since `Var(l)` and `Var_final` come from the same chains at two
  sweep indices, not independent samples as a textbook two-sample F-test
  assumes. Not supplied automatically: `n_chains` is not stored on the
  `two_block_sweep_history` object and must be pulled from the fitted
  object that produced it -- `nrow(fit$fixef[[k]])` (any RE component `k`)
  for the main-stage history, `fit$pilot_chisq$n_pilot` for the pilot-stage
  history. Omitting `n_chains` (the default) skips the band entirely, so
  existing calls are unaffected.

* **Enhancement: `plot_var_convergence()` now varies point *shape*
  (in addition to colour) across series/eigenvalue traces**, in both the
  `"base"` engine (`pch`) and the `"ggplot"` engine (`shape`, via a shared
  `scale_shape_manual()`), recycling a 12-shape palette
  (filled/open circle, triangle, square, diamond, plus, cross, asterisk).
  Colour alone (`grDevices::hcl.colors()`) can have similar hues for
  adjacent series, especially with many Block~2 fixed effects or whitened
  eigen-components; shape gives a second, easier-to-distinguish channel.

* **Bug fix: `plot_var_convergence(..., engine = "base")`'s legend
  was never displayed correctly**, in two separate ways:
    1. It was originally positioned via `par("usr")` data-range arithmetic
       (a fixed fraction of the plot's own y-range below `usr[3]`), which
       for any y-range larger than a couple of units lands well outside the
       figure/device area and is silently clipped -- i.e. the legend box
       itself was computed but rendered off-screen every time.
    2. After moving the legend into its own `graphics::layout()` panel below
       the main plot, it used a hardcoded `ncol = min(n_series, 6)` --
       `graphics::legend()` does not auto-wrap based on available width, so
       long `"re_component | covariate"` labels (e.g. `"distracted_a1 |
       free_reduced_lunch"`) didn't fit six-across on an ordinary device,
       and the overflow columns ended up positioned outside the plot area
       (e.g. "Ex_12"'s 7 Block~2 fixed effects).

  The legend now lives in its own layout panel sized to the number of rows
  actually needed, and its column count is chosen by measuring the widest
  label's rendered width against the current device width
  (`graphics::strwidth()`/`par("din")`), so it is robust to the data's
  y-range, the device size, and the length of the coefficient labels being
  plotted, not just the count of series.

* **New: `plot_var_convergence()`**, a combined (single-chart)
  `Var(l) / Var_ref` convergence plot for a `two_block_sweep_history` object's
  Block~2 fixed effects, implementing Claim~3 of the two-block Gibbs
  ergodicity reference (`inst/BLOCK_GIBBS_ERGODICITY.md`): after `l` inner
  sweeps, the cross-chain covariance \eqn{\Sigma^{(l)}_{11}} satisfies
  \eqn{\Sigma_{11}^{-1/2} \Sigma^{(l)}_{11} \Sigma_{11}^{-1/2} = I - A^{2l}}
  for a convergent matrix `A`, so `Var(l) / Var_ref` ratios are bounded above
  by 1 and increase toward it as `l` grows. Two modes:
    - Default (`whitened = FALSE`): one named-coefficient trace per
      `coef_focus` entry.
    - `whitened = TRUE`: whitens each sweep's full cross-chain covariance by
      the reference covariance and plots its (basis-invariant) eigenvalues,
      labeled `var1, var2, ...` instead of named coefficients, since the
      whitening/eigen-rotation mixes the original coefficients together.

  The reference variance/covariance is resolved automatically, not as an
  opt-in toggle: when `design` and `measurement_prior_list` are both
  supplied, the *exact* reference covariance is computed via the new
  `lmerb_posterior_covariance()` (valid only when dispersion and the
  random-effect variance-covariance are both fixed, not sampled); otherwise
  it falls back to the empirical last-sweep cross-chain covariance
  (`"Var_final"`) -- the only option for estimated-dispersion or
  estimated-vcov models.

  Two supporting pieces:
    - New exported `lmerb_posterior_covariance(design, measurement_prior_list)`:
      companion to `lmerb_posterior_mean()`, returning the exact posterior
      covariance of the stacked Block~2 hyperparameter vector (\eqn{M^{-1}},
      where `M` is the same Schur-complement posterior precision
      `lmerb_posterior_mean()` already solves against).
    - `two_block_sweep_history` objects built by the sweeps-outer/chains-inner
      engines (the `_run_with_pilot()` family: `rLMMNormal_reg_estimated_vcov()`,
      `rLMMindepNormalGamma_reg_known_vcov()`,
      `rLMMindepNormalGamma_reg_estimated_vcov()`, and -- since the
      `rGLMM_sweep()` engine swap below -- `rLMMNormal_reg_known_vcov()`'s
      `sim_method = "TWO_BLOCK_GIBBS"` route too) now also carry
      `cov_by_sweep` (the full per-sweep cross-chain covariance matrix of the
      stacked fixed effects, not just the marginal `sd` already in `table`)
      and `coef_index` (its row/column stacking order) -- purely additive
      fields, required for `whitened = TRUE`.

* **Bug fix: `Prior_Setup_GLMM(..., dispformula = ~<group>)`'s internal
  `glmmTMB` reference fit (`$fit_ref`/`$dispersion_fit`) printed its entire
  training data inline whenever a caller ran `print(summary(fit_ref))` or
  `summary(fit_ref)$call`.** Root cause: `.lmebayes_fit_glmmtmb_reference()`
  builds its `glmmTMB::glmmTMB()` call via `do.call()` with `data`/`family`
  bound to already-realized R objects (not symbols referencing a variable in
  some calling frame); `glmmTMB()`'s own `match.call()`-based `$call` then
  embeds those literal objects, and `print.summary.glmmTMB()`'s `"Data:
  "`/family header lines deparse `fit$call$data`/`fit$call$family` verbatim
  -- for a literal data frame, this dumps the whole data set as text. The
  returned `glmmTMB` fit's `$call$data`/`$call$family` are now overwritten
  with placeholder symbols (`quote(data)`/`quote(family)`) right after
  fitting; this only changes what is *displayed* for that fit's call --
  `fixef()`, `vcov()`, `predict()`, `sigma2_group`, and every other
  downstream calibration quantity are computed from the fit's actual data
  (not from `$call`) and are unaffected.

* **Bug fix: `Prior_Setup_GLMM()`'s per-RE-component `Sigma_fixef`
  (Block~2 hyperparameter prior covariance) could be numerically
  non-symmetric to floating-point precision**, causing
  `pfamily_list()`/`dGamma_list()` to fail downstream with `"matrix Sigma
  must be symmetric"` from `dNormal()`/`dIndependent_Normal_Gamma()`'s
  `isSymmetric()` validation. Root cause: `Sigma_fixef` is a column/row
  submatrix of `vcov(fit_ref)` (an `lme4`/`glmmTMB` reference fit), and
  `glmmTMB`'s Hessian-based `vcov()` is not always exactly symmetric (seen at
  the `1e-12`-ish absolute scale on a `big_word_club` `dispformula =
  ~school_id` fit with `max(abs(V - t(V)))` \eqn{\approx} \code{1.3e-12}).
  `Sigma_fixef` is now explicitly symmetrized (`(Sigma_fixef +
  t(Sigma_fixef)) / 2`) right after construction; this only removes
  floating-point noise and does not change any calibration math.

* **Bug fix: `print.Prior_Setup_GLMM()` errored on R \eqn{\ge} 4.6 when
  printing per-group Block~1 dispersion calibration** (i.e. after
  `Prior_Setup_GLMM(..., dispformula = ~<group_name>)`). It called
  `round()` on a data frame with a non-numeric `group` column, which R 4.6's
  stricter `Math.data.frame` group generic now rejects
  (`"non-numeric-alike variable(s) in data frame: group"`). Now rounds only
  the numeric columns before printing.

* **Bug fix: `rLMMNormal_reg_known_vcov()`'s exact-iid engine
  (`sim_method = "DEFAULT"`, i.e. `rLMMNormal_reg_known_vcov_iid()` /
  `rLMMNormal_joint_iid()`) returned the wrong `ranef.mode`.** It was set
  from `rLMMNormal_joint_iid()`'s `b_last` field -- the random effects
  \eqn{\beta_j} from the *last stored draw only* -- instead of a stable
  point estimate. This made `ranef.mode` from the iid route inconsistent
  with the meaning of `ranef.mode` everywhere else (the two-block Gibbs
  routes' ICM mode, which -- for these fully Gaussian-conjugate models --
  coincides with the exact posterior mean of \eqn{\beta_j}, verified against
  `mean(fit$coefficients)`); comparisons against `lme4::coef()` or against
  the mean of `fit$coefficients` could be off by many posterior standard
  deviations. `rLMMNormal_joint_iid()` now also returns `b_mean` -- the
  exact posterior mean of \eqn{\beta_j} (the group-\eqn{b_j} analogue of its
  existing `fixef_mean`), computed directly from `fixef_mean` with no Monte
  Carlo noise -- and `rLMMNormal_reg_known_vcov_iid()`'s `ranef.mode` is now
  populated from `b_mean` instead of `b_last`. `b_last` is unchanged and
  still available for callers that want the final draw specifically.
  `rLMMNormal_reg_known_vcov_two_bg()` and every other route were unaffected
  (their `ranef.mode` was already the ICM mode, not a raw draw).

* **Renamed the `x` argument to `D` and `x_hyper` to `W` on the 11
  `rLMM_reg`/`rGLMM_reg`-family matrix-level exports.** `rLMMNormal_reg()`,
  `rLMMNormal_reg_known_vcov()`, `rLMMNormal_reg_known_vcov_iid()`,
  `rLMMNormal_reg_known_vcov_two_bg()`, `rLMMNormal_reg_estimated_vcov()`,
  `rLMMindepNormalGamma_reg()`, `rLMMindepNormalGamma_reg_known_vcov()`,
  `rLMMindepNormalGamma_reg_estimated_vcov()`, `rGLMM_reg()`,
  `rGLMM_reg_known_vcov()`, and `rGLMM_reg_estimated_vcov()` now take `D`
  instead of `x` (the level-1, `l2 x p_re` random-effect design matrix) and
  `W` instead of `x_hyper` (the named list of group-level hyper-design
  matrices), matching the \eqn{D}/\eqn{\mathcal{W}} notation documented in
  `inst/notation.md` and in `?rLMM_reg`'s "Model and notation" section. This
  is purely a rename: everything previously documented for `x`/`x_hyper`
  (unique/non-empty `colnames(x)` requirement, `names(x_hyper)` matching
  `colnames(x)`, etc.) applies unchanged to `D`/`W`. `matrix_args_lmm()`'s
  and the internal `.lmebayes_matrix_args_glmm()`'s returned argument lists
  also now use `D`/`W` elements instead of `x`/`x_hyper` (matching the
  renamed formals on the routed exports they target via `do.call()`). This
  is a breaking change for any direct caller of these 11 exports passing
  `x =`/`x_hyper =` by keyword; `lmerb()`/`glmerb()`/`rlmerb()`/`rglmerb()`
  callers are unaffected. `two_block_rNormal_reg()`, `rLMMNormal_joint_iid()`,
  `two_block_rate()`, `two_block_rate_from_pfamily_list()`,
  `two_block_mode_weights()`, and `build_mu_all()` are unrelated and keep
  their existing `x`/`x_hyper` formals.

* **Renamed the grouping-factor argument `block` to `group` on all 13
  matrix-level LMM/GLMM exports.** `rLMMNormal_reg()`,
  `rLMMNormal_reg_known_vcov()`, `rLMMNormal_reg_known_vcov_iid()`,
  `rLMMNormal_reg_known_vcov_two_bg()`, `rLMMNormal_reg_estimated_vcov()`,
  `rLMMindepNormalGamma_reg()`, `rLMMindepNormalGamma_reg_known_vcov()`,
  `rLMMindepNormalGamma_reg_estimated_vcov()`, `rGLMM_reg()`,
  `rGLMM_reg_known_vcov()`, `rGLMM_reg_estimated_vcov()`,
  `two_block_rNormal_reg()`, and `rLMMNormal_joint_iid()` now take `group`
  instead of `block` as their grouping-factor formal (still required to be a
  `factor`, as above). This is purely a rename: everything below describing
  `block`'s behavior (factor requirement, `group_levels`/`group_name`
  derivation, etc.) applies unchanged to `group`; in particular
  `attr(group, "group_name")` replaces `attr(block, "group_name")` as the
  attribute checked when the caller does not pass `group` as a bare
  variable. `matrix_args_lmm()`'s returned argument list also now uses a
  `group` element instead of `block` (matching the renamed formal on the
  routed export it targets via `do.call()`). The Gibbs "two-block"
  terminology (`Block~1`/`Block~2`, `two_block_rNormal_reg()`,
  `two_block_rate()`, etc.), the generic `rNormal_reg_group()`/
  `rNormalGLM_reg_group()` block-partition family, and the compiled Rcpp/C++
  boundary are unrelated and unaffected by this rename. This is a breaking
  change for any direct caller of these 13 matrix-level exports passing
  `block =` by keyword; `lmerb()`/`glmerb()`/`rlmerb()`/`rglmerb()` callers
  are unaffected.

* **Removed the `re_coef_names` and `group_levels` arguments from all 13
  matrix-level LMM/GLMM exports; `block` must now be a factor; fixed
  `group_name` auto-derivation.** `re_coef_names` and `group_levels` were
  always derivable from `colnames(x)`/`levels(block)`, but were accepted as
  separate, never-cross-checked arguments (with a silent "synthesize RE1,
  RE2, ..." fallback when `x` had no column names). `rLMMNormal_reg()`,
  `rLMMNormal_reg_known_vcov()`, `rLMMNormal_reg_known_vcov_iid()`,
  `rLMMNormal_reg_known_vcov_two_bg()`, `rLMMNormal_reg_estimated_vcov()`,
  `rLMMindepNormalGamma_reg()`, `rLMMindepNormalGamma_reg_known_vcov()`,
  `rLMMindepNormalGamma_reg_estimated_vcov()`, `rGLMM_reg()`,
  `rGLMM_reg_known_vcov()`, `rGLMM_reg_estimated_vcov()`,
  `two_block_rNormal_reg()`, and `rLMMNormal_joint_iid()` all drop both
  formals: `x` must now have unique, non-empty `colnames(x)`, and `block`
  must now be a `factor` (there is no `group_levels` override -- use
  `factor(block, levels = full_superset)` to control level order or supply
  a level superset not present in the observed data). This also fixes a
  latent bug where `x_hyper` was silently reordered (risking `NULL`
  entries) only when `names(x_hyper)` did **not** match `colnames(x)`; it
  now errors on a mismatch and reorders only when the name sets agree.
  Separately, `group_name` auto-derivation was fixed: it previously lived
  one call-frame too deep and always resolved to the literal text
  `"group"`; it is now captured via `substitute()` in each export's own
  frame, correctly resolving to the caller's actual `block` variable name
  (erroring, instead of guessing, when `block` is not a plain variable,
  e.g. `block = df$school_id`). This is a breaking change for any direct
  caller of these 13 matrix-level exports; `lmerb()`/`glmerb()`/`rlmerb()`/
  `rglmerb()` callers are unaffected (`matrix_args_lmm()`/
  `.lmebayes_matrix_args_glmm()` already sourced these from a `design`
  object with a factor `groups` and named `Z`/`re_coef_names`, and no
  longer forward `re_coef_names`/`group_levels` at all).

* **Removed the `group_name` argument from all 13 matrix-level LMM/GLMM
  exports.** `group_name` cannot be derived from `block`'s *value* (R
  variable names are not part of an object's data), so it is now resolved
  from, in order: (1) `attr(block, "group_name")`, if set; (2)
  `substitute(block)` in the caller's own frame, when `block` is passed as
  a bare variable (e.g. `block = school_id`). Direct callers that pass
  `block` as a non-symbol expression (e.g. `block = df$school_id`) must now
  attach the name themselves via
  `attr(block, "group_name") <- "school_id"` beforehand, instead of passing
  a separate `group_name` argument; calling with an explicit `group_name =`
  argument now errors with "unused argument". `lmerb()`/`glmerb()`/
  `rlmerb()`/`rglmerb()` callers are unaffected (`matrix_args_lmm()`/
  `.lmebayes_matrix_args_glmm()` now attach the attribute internally
  instead of forwarding `group_name`).

* **Removed the required `P`/`Sigma` field from `prior_list`/
  `prior_list_block1` for `rGLMM_reg()`, `rGLMM_reg_known_vcov()`,
  `rGLMM_reg_estimated_vcov()`, `two_block_rNormal_reg()`, and
  `rLMMNormal_joint_iid()`.** Extends the `rLMMNormal_reg`-family `P`
  removal (below) to these five exports: the Block~1 random-effect prior
  precision was always exactly `diag(tau2_k)` from `pfamily_list` (the same
  plug-in `priors_from_pfamily_list()` already computes), but was required
  as a `prior_list$P`/`prior_list$Sigma` field, never cross-checked against
  `pfamily_list`. These five exports now derive it internally via
  `.rLMM_P_from_pfamily_list()` right after validating `pfamily_list`/
  `re_names`, and `stop()` if a caller supplies `prior_list$P`/`$Sigma` (or
  `prior_list_block1$P`/`$Sigma`) themselves. This is a breaking change for
  any direct caller of these five exports; `lmerb()`/`glmerb()`/`rlmerb()`/
  `rglmerb()` callers are unaffected -- `.lmebayes_block1_prior_list()` no
  longer builds a `P` either (dropping the informational, and since the
  earlier `matrix_args_lmm()` fix, stale `P` previously shown in
  `rlmerb()`/`rglmerb()`'s returned `$Prior$block1_prior` field).
  `two_block_rate()`/`two_block_rate_from_pfamily_list()` (diagnostic rate
  calibration, not a regression engine) and `rGLMM_sweep()` are unaffected
  and still require `P`/`Sigma` directly.

* **Removed the `P` argument from all eight `rLMMNormal_reg`/
  `rLMMindepNormalGamma_reg` exports.** `P` (the Block~2 random-effect prior
  precision matrix) was always mechanically derivable from `pfamily_list`
  (one \eqn{\tau^2_k} plug-in per component -- fixed `dispersion` for
  `dNormal()`, prior mean `rate/(shape - 1)` for
  `dIndependent_Normal_Gamma()`), but was accepted as a separate,
  never-cross-checked argument, so a caller could silently pass a `P`
  inconsistent with `pfamily_list`. `rLMMNormal_reg()`,
  `rLMMNormal_reg_known_vcov()`, `rLMMNormal_reg_known_vcov_iid()`,
  `rLMMNormal_reg_known_vcov_two_bg()`, `rLMMNormal_reg_estimated_vcov()`,
  `rLMMindepNormalGamma_reg()`, `rLMMindepNormalGamma_reg_known_vcov()`, and
  `rLMMindepNormalGamma_reg_estimated_vcov()` all drop the `P` formal; each
  now derives it internally via a new `.rLMM_P_from_pfamily_list()` helper
  right after validating `pfamily_list`. This is a breaking change for any
  direct caller of these matrix-level exports (`lmerb()`/`glmerb()`/
  `rlmerb()`/`rglmerb()` callers are unaffected -- `matrix_args_lmm()` no
  longer builds or forwards `P` either). `rLMMNormal_joint_iid()` and the
  GLMM Block~1 `prior_list_block1$P` shape are unrelated and unchanged.

* **Exported (and renamed) two internal helpers used by `lmerb()`/`glmerb()`:
  `priors_from_pfamily_list()` (was `.lmebayes_priors_from_pfamily_list()`)
  and `matrix_args_lmm()` (was `.lmebayes_matrix_args_lmm()`).** Both were
  previously only reachable from **lmebayes** via `:::`; they are now
  ordinary exports (marked `@keywords internal`, so they stay out of the
  default help index) with `@param`/`@return` documentation. Their
  argument lists and return shapes are unchanged and are still considered
  refactor candidates -- both currently bundle several responsibilities
  (dispersion resolution, `pfamily_list` validation, route-specific argument
  branching) that will likely be simplified/split in a follow-up; treat
  the current signatures as provisional. All in-package and **lmebayes**
  call sites were updated to the new names.

* **New `sim_method` argument: exact iid sampling for the fixed-dispersion /
  known-variance-components route.** `rlmerb()`, `rglmerb()` (for
  `family = gaussian()`), and `rLMMNormal_reg_known_vcov()`/`rLMMNormal_reg()`
  gain a `sim_method` argument: `"DEFAULT"` or `"TWO_BLOCK_GIBBS"`. This only
  changes behavior on the `lmm_fixed_known` route -- fixed (scalar or
  per-group) observation dispersion **and** every Block~2 component
  `dNormal()` (known/fixed variance components) -- where the joint posterior
  over `(gamma, b_1, ..., b_J)` is exactly multivariate normal (see
  `inst/README_KNOWN_VCOV_GAUSSIAN.md`). `sim_method = "DEFAULT"` now draws
  directly, iid, from that closed-form Gaussian via a new
  `rLMMNormal_reg_known_vcov_iid()`/`rLMMNormal_joint_iid()` engine -- no
  Gibbs sweeps, no burn-in, no residual autocorrelation between draws.
  `sim_method = "TWO_BLOCK_GIBBS"` keeps the previous two-block Gibbs engine,
  now factored out unchanged as `rLMMNormal_reg_known_vcov_two_bg()`.
  `rLMMNormal_reg_known_vcov()` itself becomes a thin dispatcher between the
  two. Every other route (`dGamma()`/`dIndependent_Normal_Gamma` Block~2
  components, or estimated variance components via
  `rLMMNormal_reg_estimated_vcov()`) only has a two-block Gibbs engine, so
  `sim_method` is a no-op there (both values behave identically). All
  matrix-level exports and the sampler output (`sim_method_used`,
  reflecting whichever engine actually ran) gain this field; no existing
  argument was renamed or removed.

  New building blocks behind this: `.lmerb_posterior_normal_system()` and
  `.lmerb_posterior_b_given_gamma()` (extracted from `lmerb_posterior_mean()`
  without changing its signature/return) build and back-substitute the exact
  joint Gaussian posterior system; `.lmerb_posterior_system_cholesky()`
  Cholesky-factors that system once per call (with a symmetry check/tolerance
  on the Block~2 precision matrix `M` that `lmerb()`/`glmerb()` always
  satisfy by construction, but direct matrix-level calls with an inconsistent
  `Sigma_ranef`/`pfamily_list` may not).

* **New `dispersion_ranef` shape: fixed per-group dispersion vector.**
  `.lmebayes_resolve_dispersion_ranef()` (used by `rlmerb()`, and by
  `rglmerb()` for `family = gaussian()`) now accepts a plain numeric vector
  of length > 1 for `dispersion_ranef`, resolved to a new
  `dispersion_mode = "fixed_vector"`: a named vector with one known, fixed
  \eqn{\sigma^2_j} per group (names must match `levels(design$groups)`
  exactly; reordered to `group_levels` order internally). Unlike the
  existing `"gamma_list"` mode (a per-group `dGamma()` *prior*, sampled and
  requiring full column rank for its ING envelope), a fixed vector is a
  directly user-supplied constant -- no sampling, no rank restriction, and
  no reference-fit dependency. Dispatches to the same
  `rLMMNormal_reg_known_vcov()`/`rLMMNormal_reg_estimated_vcov()` routes as
  the existing scalar `"fixed"` mode (`.rLMM_validate_fixed_dispersion_prior_list()`
  now accepts `prior_list$dispersion` as a length-`J` vector, threaded
  through a new `group_levels` parameter). The C++ engine
  (`rNormalRegBlocks()`) already supported a per-group dispersion vector;
  no `src/*.cpp` changes were needed.

  Fixed two latent scalar-dispersion assumptions in the Gaussian rate-
  calibration and ICM code paths, uncovered while verifying the new mode
  end to end (both previously silently used only the *first* group's
  dispersion for every group when handed a vector): the per-observation
  working-weight derivation in `.two_block_rate_inputs()` (Theorem~3
  convergence-rate calibration) now expands a per-group `dispersion`
  vector via the same group-to-row mapping used elsewhere in that
  function; and `lmerb_posterior_mean()`/`glmerb_posterior_mode()` now
  index `measurement_prior_list$dispersion_ranef` per group inside their
  per-group loops instead of dividing/passing the whole vector at once.

* **Stage 3a/3b — dead iid C++ removal (deduplication):** Following the
  Stage 1a/1b/1c R-level removals, audited `lmebayesCore`'s C++ for iid
  routines that are now unreachable from R. Confirmed `rNormalGammaReg()`,
  `rGammaGaussian()`, and `rGammaGamma()` (in `src/rNormalGammaReg.cpp`,
  `src/rGammaGaussian.cpp`, `src/rGammaGamma.cpp`) were dead: their only
  callers were `rNormalGamma_reg()`, `rGamma_reg()`, and
  `rGamma_Conjugate_reg()` in `R/simfunction.R`, which Stage 1b already
  removed (those live on in `glmbayesCore`). Deleted the three `.cpp` files,
  their `_cpp_export()` wrappers in `src/export_wrappers.cpp`, the matching
  declarations in `src/simfuncs.h`, and the internal `.rNormalGammaReg_cpp()`
  / `.rGammaGaussian_cpp()` / `.rGammaGamma_cpp()` `.Call()` shims in
  `R/rcpp_wrappers.R`; regenerated `RcppExports.R`/`RcppExports.cpp` via
  `Rcpp::compileAttributes()`. No NAMESPACE or public-API change (all were
  internal, `@noRd`). The remaining mixed-model C++ (`twoBlockGibbs.cpp`,
  `block_rIndepNormalGammaReg.cpp`, `rNormalRegBlocks.cpp`,
  `rNormalGLMBlocks.cpp`, etc.) still directly links against the rest of the
  iid envelope/sampler engine under `src/`, so most of that engine remains
  compiled into this package; pruning it requires a `glmbayesCore`
  C-callable bridge and is deferred to Stage 3c–3f pending `glmbayesCore`'s
  CRAN review.

* **Stage 1a/1b/1c — remove duplicate iid material; delegate to `glmbayesCore`
  (deduplication):** Following Stage 0's dependency wiring, removed 21 R files
  duplicating `glmbayesCore` functionality (~47 exports total: truncated-dist
  and envelope C++ callbacks, the `simfunction`/`multi_r*` iid sampler stack,
  and the top-level `Prior_Setup`/`pfamily`/`rglmb`/`rlmb` prior API), along
  with their 32 `man/` pages and 28 `inst/examples/` files. **Hard break, no
  aliases:** callers that used `lmebayesCore::Prior_Setup()`,
  `lmebayesCore::rglmb()`, `lmebayesCore::dNormal()`, etc. must switch to
  `glmbayesCore::…` directly. Remaining internal call sites in
  `lmebayesCore`'s own R/ (two-block engines, `Prior_Setup_GLMM()`,
  `lmebayes_posterior_icm.R`, etc.) were updated to call `glmbayesCore::…`
  explicitly. `lmebayesCore` now exports **42 symbols + 6 S3 methods**
  (mixed-model setup and two-block Gibbs sampling only).
  `rindepNormalGamma_reg_with_envelope()` -- caught by the mass file removal
  but actually `lmebayesCore`-specific, not a `glmbayesCore` duplicate -- was
  restored as its own file.

* **Stage 2 — C++ R-namespace retarget:** `src/package_ns.h`'s `GLMBAYES_R_NS`
  macro (used by `src/R_interface.h` and `src/twoBlockGibbs.cpp` to resolve
  R-level callbacks -- `EnvelopeOpt`, `EnvelopeSort`, `glmbfamfunc`,
  `rNormal_reg.wfit`, `rgamma_ct`, and Block~2's `rglmb`) now points at
  `"glmbayesCore"` instead of `"lmebayesCore"`, matching the R-level exports
  removed above. Kept as a single macro rather than splitting into separate
  iid (`glmbayesCore`) / mixed-model (`lmebayesCore`) namespace macros, since
  nothing in the remaining mixed-model C++ needs a purely `lmebayesCore`-local
  R callback today. **C++ object code itself is not yet deduplicated** -- the
  full iid envelope/sampler engine under `src/` is still compiled into this
  package alongside the mixed-model-only `.cpp` files that link against it
  (Block~1 in `twoBlockGibbs.cpp`, `block_rIndepNormalGammaReg.cpp`,
  `rNormalRegBlocks.cpp`, `rNormalGLMBlocks.cpp`); pruning that duplication is
  a future Stage 3 effort.

* **Stage 0 — `glmbayesCore` dependency wiring (deduplication prep):**
  `DESCRIPTION` now `Imports: glmbayesCore (>= 0.5.1)`. `.onLoad()` fails
  fast with an install hint if `glmbayesCore` is missing. No exports or
  duplicate sources removed yet.

* **Forked from `glmbayesCore` as the full-featured backend for `lmebayes`:**
  `lmebayesCore` is a history-preserving fork of `glmbayesCore` (created
  2026-07-15) that keeps the complete glm/envelope engine *and* the
  two-block Gibbs mixed-model stack (`model_setup()`, `Prior_Setup_GLMM()`,
  `rlmerb()`/`rglmerb()`, `rLMM_reg`/`rGLMM_reg` routes, etc.). `glmbayesCore`
  itself is being stripped down to only the glm/envelope engine that
  `glmbayes` needs; `lmebayes` now depends on `lmebayesCore` instead.
  Package identity (DESCRIPTION, `NAMESPACE`, `GLMBAYES_R_NS`,
  `.Call()`/DLL registration) was renamed accordingly; the internal C++
  `glmbayes::` namespace was left unchanged as an implementation detail.
  First time builds


* **Per-group `dGamma_list()` prior now uses the §3.3.4 marginal rate:**
  **`dGamma_list.Prior_Setup_GLMM()`** feeds each group's `dGamma()` the
  Chapter A12 **§3.3.4** marginal ING rate (`beta` integrated out) instead of
  the **§3.3.5** fixed-`beta` `rate_gamma`. This is the theoretically correct
  choice for the Block~1 ING sampler, which draws `sigma2_j` from the
  marginal law and then `b_j | sigma2_j` (`beta` is never held fixed at a
  point estimate during that draw). `rate_gamma` remains on
  `ing_prior_measurement_group` for diagnostic comparison only (printed by
  a dev-only table in **`Prior_Setup_GLMM()`** whenever `dispformula`
  requests per-group dispersion); nothing downstream consumes it. Truncation
  bounds (`disp_lower`/`disp_upper`, `blup_infl`, `R_lo`/`R_hi`) are
  unaffected -- they were already mean-matched at `sigma2_hat`, which is now
  also the new rate's exact prior mean (previously it was not, for the
  most BLUP-inflated groups). See `inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md`
  Parts I-II.

* **GLMM router (`.lmebayes_run_glmm_engine()`):** non-Gaussian **`rglmerb()`**
  dispatches through **`REG_ROUTE_TABLE`** to **`rGLMM_reg_known_vcov()`** or
  **`rGLMM_reg_estimated_vcov()`** (replacing a direct **`rGLMM_reg()`** call).
  LMM routing (**`.lmebayes_run_lmm_engine()`**) uses the same table for four
  Gaussian routes.

* **GLMM engines split (`rGLMM_reg.R`):** monolithic **`rGLMM()`** replaced by
  **`rGLMM_reg_known_vcov()`**, **`rGLMM_reg_estimated_vcov()`**, and dispatcher
  **`rGLMM_reg()`** (shared help **`?rGLMM_reg`**). Non-Gaussian models always
  run a pilot stage (unless **`n_pilot = 0L`**); the two routes differ in
  eigenvalue-bound complexity (fixed **`dNormal`** τ² vs ING **`disp_lower`**
  conservatism), not in whether a pilot runs.

* **LMM engines merged (`rLMM_reg.R`):** **`rLMMNormal_reg.R`** and
  **`rLMMIngNormal_reg.R`** are one module with shared help **`?rLMM_reg`**
  (aliases for six exports). There is no standalone **`rLMM()`** export —
  matrix Gaussian LMMs use **`rLMMNormal_reg_*`** / **`rLMMindepNormalGamma_reg_*`**
  routes; formula GLMMs use **`rGLMM_reg`** via **`rglmerb()`**.

* **`rlmerb()` four-route routing:** **`.lmebayes_run_lmm_engine()`** dispatches
  to **`rLMMNormal_reg_known_vcov()`**, **`rLMMNormal_reg_estimated_vcov()`**,
  **`rLMMindepNormalGamma_reg_known_vcov()`**, or
  **`rLMMindepNormalGamma_reg_estimated_vcov()`** from fixed vs ING Block~2 and
  fixed vs dGamma σ². Legacy **`rLMMindepNormalGamma_reg()`** (outer σ² loop)
  remains exported but is not the default **`rlmerb()`** path.

* **ICM at fixed variance components:** **`lmerb_posterior_mean()`** /
  **`glmerb_posterior_mode()`** iterate Block~1 / Block~2 hyperparameters at
  **fixed** \(\tau^2_k\) and \(\sigma^2\) plug-ins
  (**`.two_block_tau2_plug_in_from_pfamily()`**, `rate/(shape−1)`). Removed
  joint posterior-mode τ² iteration (**`two_block_joint_posterior_mode()`** stack).

* **Scale-invariant ICM stopping rule:** **`glmerb_posterior_mode()`** now
  measures the Block~2 `delta` used for the `tol`/`converged` check as a
  per-component Mahalanobis distance \(\sqrt{(\gamma_k^{new}-\gamma_k)^\top
  P_{\gamma_k}^{\mathrm{post}} (\gamma_k^{new}-\gamma_k)}\) in each RE
  component's own posterior-precision metric, maximized over components
  \(k\) — not the previous raw \(\ell_\infty\) change in `fixef`. The old
  criterion depended on the arbitrary units of each hyper-covariate, so
  rescaling or whitening an `X_hyper[[k]]` column changed convergence
  behavior even though the fitted model and posterior were unchanged; the
  new one does not. The per-block posterior precision (`post_P_list`) is now
  also hoisted out of the ICM loop since it does not depend on the iteration
  state. (Superseded for **`lmerb_posterior_mean()`** by the exact closed-form
  solve below, which has no stopping rule at all.)

* **`lmerb_posterior_mean()` is now an exact closed-form solve, not ICM:**
  the Gaussian Block~1/Block~2 target is exactly jointly Gaussian, and
  Block~1's conditional mean per group is affine in the shared hyperparameter
  vector \(\gamma\) with no direct coupling between groups. Substituting that
  affine relationship into the Block~2 update eliminates every group's random
  effect algebraically (a Schur-complement/Henderson-mixed-model-equations
  elimination), leaving one small linear system in \(\gamma\) alone (dimension
  = total hyperparameter count, independent of the number of groups \(J\)).
  Solving it once gives the exact joint mean — no alternating iteration,
  `tol`, `maxit`, or non-convergence warning is possible for this model, and
  it costs \(O(J)\) (never a \(J \times J\) or \(J p_{re}\)-dimensional
  matrix), so it scales to large numbers of groups. `tol`/`maxit` remain
  accepted (for interface parity with **`glmerb_posterior_mode()`**, which is
  unchanged and still iterates for non-Gaussian families) but are unused;
  the return always has `converged = TRUE`, `iterations = 1L`, `delta = 0`.
  This also restores the exact `D0 = 0` (start at the true posterior mean)
  assumption that **`two_block_tv_bound()`**/**`two_block_l_for_tv()`**'s
  sweep-count guarantee relies on, which a non-converged ICM start could
  silently violate.

* **`two_block_l_for_tv()` no longer errors when the search exceeds `l_max`:**
  it now issues a single `warning()` (when `warn = TRUE`) and returns a
  practical uncertified fallback capped at **200 inner sweeps**
  (`l <= 199`), not `l_max = 1e6`.  The old `l_max` return caused integer
  overflow in pilot cost optimization (`n_pilot * m_convergence_pilot`)
  and invalid `n_pilot` values.  Internal repeated calls (pilot cost
  search) pass `warn = FALSE` to avoid warning spam.  Mode-gap pilot sweep
  calibration and inner-sweep counts are likewise capped at 200 via
  `.two_block_cap_inner_sweeps()` / `.two_block_m_pilot_from_gap()`.
  `l_max`/`m_min` calibration is inherently a best-effort setup step
  (choosing a burn-in sweep count), not part of the returned draws, so a
  near-degenerate `rate$lambda_star` (close to 1) should not abort the
  whole **`lmerb()`**/**`rlmerb()`**/**`glmerb()`**/**`rglmerb()`**
  call; all internal calibration call sites inherit this automatically.

* **`dGamma_list()` gains `disp_center = c("sigma2_hat", "dispersion2")`:**
  the default (`"sigma2_hat"`) reproduces the existing per-group truncation
  window unchanged (mean-matched at `sigma2_hat_j`, upper tail widened by
  `disp_upper_anchor`). The new opt-in `"dispersion2"` mean-matches both
  bounds symmetrically at an `EnvelopeCentering()`-style dispersion estimate
  that integrates over the random effect's own posterior uncertainty
  (new `n_rss_iter` argument, default `10L`, controls its fixed-point
  iteration count) instead of BLUP-inflating `sigma2_hat_j`; `disp_upper_anchor`
  is ignored in this mode. Tends to produce narrower, better-centered upper
  tails for groups with large BLUP/OLS RSS inflation. New helper
  **`.lmebayes_group_dispersion2_envelope_centering()`**
  (`R/mixed_rmerb_helpers.R`); `window_diagnostics` gains a `dispersion2`
  column. See `inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md` Part III.

* **Rate helper rename:** **`two_block_rate_v2()`** removed; use
  **`two_block_rate_from_pfamily_list()`** (`R/two_block_ergodicity.R`)
  for the `pfamily_list` adapter around **`two_block_rate()`**.

* **Two-block ergodicity consolidation:** rate, TV-bound, and mode-weight
  helpers merged into **`R/two_block_ergodicity.R`**. **`two_block_mode_weights()`**
  is no longer exported (still used internally by **`rGLMM()`** for non-Gaussian
  rate calibration).

* **C++ R callbacks via registered namespace:** envelope and simulation C++
  now resolve **`EnvelopeOpt()`**, **`EnvelopeSort()`**, **`glmbfamfunc()`**,
  **`rNormal_reg.wfit()`**, and **`rgamma_ct()`** from the **`glmbayesCore`**
  namespace (`R_interface.h` / `GLMBAYES_R_NS`), so downstream packages
  (e.g. **lmebayes**) no longer need to re-export them for search-path lookup.

* **LMM engine split:** **`rLMMNormal_reg()`** samples with fixed observation
  dispersion (\code{prior_list = list(dispersion = sigma2)}); **`P`** is a
  separate argument. **`rLMMindepNormalGamma_reg()`** implements an outer
  two-block Gibbs sampler (dispersion via **`rGamma_reg()`**, fixed effects and
  random effects via full **`rLMMNormal_reg()`** runs).

* **`rGLMM()` pilot defaults:** non-Gaussian models now run pilot + main by
  default. New argument **`gap_tol`** (default `0.0196`) derives **`n_pilot`**
  when **`n_pilot = NULL`**; **`n_pilot = 0L`** or **`gap_tol = NULL`** skips
  the pilot. Gaussian models never run a pilot. **`tv_tol`** now defaults to
  **`0.01`**. Helper **`.two_block_resolve_n_pilot()`** centralises the policy.

* **Matrix LMM / GLMM ICM:** when `start = NULL` (default), **`rLMM_reg`**
  routes and **`rGLMM()`** compute Block~2 starts via
  **`lmerb_posterior_mean()`** / **`glmerb_posterior_mode()`** at **fixed**
  variance-component plug-ins, using **`.two_block_measurement_prior_list()`**
  (and **`.two_block_tau2_plug_in_from_pfamily()`** for τ²). Outputs include
  **`ranef.mode`** and **`icm_info`**. Non-Gaussian **`rGLMM()`** still requires
  **`b_start`** when **`start`** is user-supplied.

* **Matrix LMM replicate chains:** Gaussian LMM sampling is exported as six
  **`rLMM_reg`** engines (four direct **`rlmerb()`** routes plus two
  dispatchers). Formula-level fitting remains in **lmebayes**
  (`rlmerb()` / `lmerb()`).

* **Restored `rGLMM()`:** matrix-level GLMM replicate-chain orchestration
  (TV calibration, pilot chi-squared, post-pilot eigenvalue upper bound,
  main-stage sampling via **`rGLMM_sweep`**) is exported again
  as **`rGLMM()`**. Replaces the earlier C++-staged implementation; returns
  the `fixef.*` namespace. Formula-level fitting remains in **lmebayes**
  (`rglmerb()` / `glmerb()`).

* **Candidate counts surfaced by the two-block v2 sampler:**
  **`two_block_rNormal_reg_v2()`** now returns `iters_fixef_draws`, an
  `n x p_re` matrix of the total number of Block 2 candidates generated per
  stored draw, summed over the `m_convergence` inner sweeps.
  `dIndependent_Normal_Gamma` components count the envelope accept-reject
  candidates until acceptance (the `iters_out` already produced by
  `rIndepNormalGammaReg`, previously discarded by the Gibbs loop); `dNormal`
  components count exactly one conjugate draw per sweep, so their columns
  equal `m_convergence`.  Dividing by `m_convergence` gives the average
  number of candidates per accepted draw (roughly the reciprocal envelope
  acceptance rate), matching the `iters` semantics of `rglmb`-style
  samplers.  Reading the counts consumes no RNG, so draws are
  bitwise-identical to the previous version under the same seed.

* **Prior-vs-data guard for `dIndependent_Normal_Gamma` sampling:**
  **`rindepNormalGamma_reg()`** now rejects calls where the Gamma (precision)
  part of the prior carries more effective prior observations than the data
  supply: inverting the `Prior_Setup()` calibration
  `shape = (n_prior + 1 + p)/2`, sampling requires
  `n_prior <= n_w = sum(weights)` (equivalently a prior weight
  `pwt <= 0.5`). Rationale: the dispersion envelope caps its log-tilt at
  `n_w/2` - the *data* contribution to the posterior Gamma shape (Remark
  4.1.3 of the ING vignette) - a strengthening of the validity condition
  `lm_log2 < shape2` that presumes a likelihood-dominated regime.
  Prior-dominated calls could previously bind that cap on every envelope
  build (console `UB3A mean slope` warnings) and silently degrade the
  envelope. Note that `n_prior` here is the effective sample size of the
  Gamma component specifically; under the `Prior_Setup()` calibration the
  Gamma and coefficient parts share a common `n_prior`, so the two are not
  fully independent.

* **Same guard in the two-block v2 sampler:**
  **`two_block_rNormal_reg_v2()`** enforces `n_prior <= J` per
  `dIndependent_Normal_Gamma` component (with `J = length(group_levels)`,
  the Block 2 hyper-regression observation count and `q_k = length(mu)`:
  `2*shape - 1 - q_k <= J`, i.e. `pwt_disp <= 0.5`).  Calibration-only
  paths (`two_block_rate_v2()`) are exempt since they use the `disp_lower`
  plug-in without sampling.

* **`pfamily_list()` generic:** New S3 generic for building a named list of
  pfamily objects from a prior-specification container.  Downstream
  packages provide methods (e.g. `lmebayes` for `Prior_Setup_GLMM()`
  objects, mapping each random-effect component to `dNormal()` or
  `dIndependent_Normal_Gamma()`).

* **Convergence rate for the two-block sampler:** New **`two_block_rate()`**
  computes the eigenvalues of
  `A = P11^{-1/2} P12 P22^{-1} P21 P11^{-1/2}` (Nygren 2020, Remark 8) for
  the joint Gaussian posterior targeted by **`two_block_rNormal_reg()`**,
  without ever forming the `J*p_re x J*p_re` Block 1 precision: the cross
  moment is accumulated per group with `p_re x p_re` solves followed by a
  single `q x q` symmetric eigendecomposition. The maximal eigenvalue
  `lambda*` is the geometric TV contraction rate of the sampler;
  `m_for_tol(tol)` returns the implied number of inner Gibbs sweeps. For
  non-Gaussian families explicit IRLS-style `weights` give a local-Gaussian
  heuristic. Validated against a dense brute-force construction of the joint
  precision and against the observed contraction of the ICM mean recursion
  (`lmerb_posterior_mean()`), which contracts at exactly `lambda*`.

* **Likelihood precision at the posterior mode:** New
  **`two_block_mode_weights()`** evaluates per-observation likelihood
  precisions (IRLS/Fisher weights) at a supplied random-effects value -
  typically the joint posterior mode from `glmerb_posterior_mode()` - and
  assembles the per-group likelihood precision blocks `Z_j' W_j Z_j`.
  Weights are computed generically from the family object
  (`w_i = wt_i mu'(eta_i)^2 / (V(mu_i) phi)`): exact observed Hessian for
  canonical links (gaussian, poisson-log, binomial-logit), expected (Fisher)
  information otherwise - including correct probit/cloglog/Gamma-log weights
  where `glmbfamfunc()$f7` carries copy-pasted logistic weights.  The
  `weights` component feeds `two_block_rate(weights = )` directly, providing
  the local-Gaussian heuristic input for extending the TV-rate analysis to
  non-Gaussian `glmerb` models.  Validated against `f7` on its correct
  branches and against the exact Gaussian rate path.

* **Explicit TV convergence bounds:** New **`two_block_tv_bound()`**
  evaluates the total-variation bound between the `l`-step kernel and the
  target (Nygren 2020) from the `two_block_rate()` spectrum, two ways:
  `method = "theorem3"` computes the exact per-eigendirection terms
  `d_i^(l)` using the closed form `erf_n(x) = pchisq(2 x^2, n)` with
  `r_i^(l) = (1 - a_{i-1}^{2l})/(1 - a_i^{2l})`; `method = "corollary1"`
  evaluates the looser geometric envelope with explicit constants. With the
  chain started at the exact posterior mean (as `lmerb` does), the mean term
  vanishes identically (`D0 = 0` default) and only the variance-convergence
  sum remains, which decays like `lambda*^{2l}` - twice the exponent of the
  crude `(lambda*)^m` proxy. **`two_block_l_for_tv()`** inverts the bound to
  give the number of inner Gibbs sweeps required for a target tolerance, and
  `print.two_block_rate()` now tabulates proxy vs Theorem 3 vs Corollary 1
  sweeps. On the lmerb big_word_club example (`lambda* = 0.839`): TV <= 1e-3
  needs 16 sweeps (Theorem 3) / 23 (Corollary 1) vs 40 for the proxy.

* **Two-block Gibbs loop in C++:** The main loop of
  **`two_block_rNormal_reg()`** (Block 1 random-effects update, Block 2
  hyperparameter update, `m_convergence` inner steps, replicate sampling) now
  runs entirely in C++ (`two_block_rNormal_reg_cpp_export` in
  `src/twoBlockGibbs.cpp`), eliminating per-iteration R/C++ round trips. This
  is a port-only change: the R wrapper still performs input validation,
  `glmbfamfunc()` resolution, and output assembly, and the C++ driver calls
  the same per-block samplers (`rNormalGLM` envelope sampler, `rNormalReg`)
  in the same order as the previous R loop. Draws are statistically
  equivalent but not bit-reproducible against the old R loop because the C++
  rejection sampler uses its own RNG stream (compare averages over many
  draws, not individual draws).

* **Faster GLM block sampling:** **`rNormalGLM_reg_group()`** now performs block
  partitioning and prior payload assembly in C++
  (`block_rNormalGLM_cpp_export`), removing per-call R overhead in block
  Gibbs loops (e.g. Block 1 of the **lmebayes** two-block sampler). The
  sampling algorithm itself is unchanged: each block still calls the existing
  `rNormalGLM()` envelope sampler serially. Posterior modes are numerically
  identical to the previous R-prep path; individual draws follow the same
  distribution but are not bit-reproducible against the old path (compare
  means over longer runs). Present-but-`NULL` prior elements (e.g.
  `dispersion = NULL`) are treated as absent, matching R `is.null()`
  semantics.

# glmbayes 0.9.6

## Highlights

* **Multi-response `lmb()`:** **`lmb()`** now handles both univariate and
  multivariate responses with a single unified interface, mirroring the behaviour
  of R's **`lm()`**. When the response has a single column the result is an
  **`lmb`** object (unchanged from prior releases). When the formula specifies
  multiple response columns (e.g. `cbind(y1, y2) ~ x`), **`lmb()`** fits a
  separate Bayesian linear model per response column and returns a named list
  with class **`mlmb`**. For the multi-response case, **`pfamily`** must be a
  list of **`pfamily`** objects with exactly one entry per response column;
  passing a single **`pfamily`** object is an error. Summary, print, and
  coefficient methods for **`mlmb`** objects are included.

* **Conjugate GLM priors (Poisson, binomial, Gamma):** New closed-form IID
  sampling paths for intercept-only models with identity links. **`dBeta()`**
  with **`rBeta_reg()`** supports Beta–Binomial(identity) conjugate updates;
  **`dGamma(Inv_Dispersion = FALSE)`** with **`rGamma_Conjugate_reg()`**
  supports Gamma–Poisson(identity) and Gamma–Gamma(identity) rate priors.
  **`Prior_Setup()`** can calibrate conjugate hyperparameters for these
  families (weighted Poisson rate and binomial probability defaults). See
  **`?dBeta`**, **`?dGamma`**, and the Chapter 02 / Chapter 07–11 vignettes.

* **Vignette structure:** Reworked **Chapter 00** as a roadmap across five
  main parts plus technical appendices. **Chapter 02** is now a conceptual
  introduction to single-parameter conjugacy; worked examples move to
  **Chapter 02-S01** through **Chapter 02-S05** (Beta–Binomial, Normal–Normal,
  Gamma–Poisson, exposure-weighted Poisson, and related topics). A **Companion
  textbooks** section in Chapter 00 indexes optional Bayes Rules! and `LearnBayes`
  appendices tied to the main GLM chapters.

* **`opencltools` import:** Core host/runtime OpenCL discovery and diagnostics
  (`detect_*`, PATH helpers, environment checks) now live in the **`opencltools`**
  package (`Imports`, >= 0.8.0). **glmbayes** keeps package-specific entry
  points (`glmbayesCore_has_opencl()`, `diagnose_glmbayes()`) that report compile-time
  OpenCL status for this build while delegating shared GPU/runtime checks—reducing
  duplicated maintenance in **glmbayes**.

* **Bayes Rules! companion examples:** Optional vignette appendices reproduce
  book datasets and published posterior summaries using **`lmb()`**, **`glmb()`**,
  **`Prior_Setup()`**, and **`dNormal()`** (suggested package **`bayesrules`** for
  data only). Coverage includes **`bikes`** (Ch. 03), **`weather_perth`** (Ch. 08–09),
  **`equality_index`** (Ch. 10), Gamma–Poisson conjugacy (Ch. 02-S04), and a
  scope note for Gamma regression (Ch. 11). Comparison tables use **printed book
  values**, not live **`rstanarm`** fits. See **Chapter 00** § Companion textbooks.

* **`LearnBayes` examples:** **Chapter 02-S04**, Appendix A, maps the
  **`hearttransplants`** example from Albert (2009) / `LearnBayes` (exposure-weighted
  Gamma–Poisson conjugacy) to **`glmb()`** with analytic Albert posteriors for
  verification (suggested package **`LearnBayes`**).

## Other changes

* Expanded **testthat** coverage for **`dBeta()`** / binomial(identity) conjugate
  paths and related **`glmb()`** integration.

# glmbayes 0.9.5

* **Tests / CRAN:** All **OpenCL**-specific **testthat** blocks now call
  **`skip_on_cran()`** (in addition to **`skip_if_no_opencl()`**), consistent
  with existing Boston/Cleveland OpenCL tests. OpenCL coverage remains for local
  checks and source builds with OpenCL; CRAN checks avoid parallel/GPU-heavy
  tests that could trigger **CPU time vs elapsed time** NOTES.

# glmbayes 0.9.4

* **Vignettes:** A vignette that previously used the `notangle` engine now
  uses the standard R Markdown vignette machinery (`knitr` /
  `rmarkdown::html_vignette`), so builds align with CRAN expectations and
  vignette index ordering should be consistent with the rest of the package.

* **OpenCL sources (`inst/cl`):** Removed unused or superseded material,
  consolidated kernels and library fragments, and aligned `.cl` layout and
  dependency tagging with the conventions used in 'openclport' and
  'nmathopencl' (prelude, shims, `nmath/` stems, family kernels under
  `src/`). See `inst/cl/README.md` for how the assembled program is stitched.

* **OpenCL program assembly:** Reworked loading so the full OpenCL program is
  built from explicit fragments (global header, `nmath` closure, family/link
  kernels) rather than ad hoc concatenation—clearer ownership of what enters
  GPU compilation and easier parity with CPU paths.

* **Tests:** Added and expanded **testthat** coverage aimed at OpenCL code
  paths (including binomial examples that exercise GPU envelope evaluation),
  complementing existing Cleveland-style checks.

* **Bug fix — binomial OpenCL:** Binomial `f2_f3` OpenCL kernels now evaluate
  the data log-likelihood with the same **proportion × trial-count**
  semantics as **`dbinom_glmb`** on the CPU (`round` successes and trials,
  clamped probability). This fixes envelope / PLSD failures for aggregated
  binomial data (e.g. `cbind(successes, failures)` / `MASS::menarche`) where
  the previous kernels treated **`y`** like a raw success count.

# glmbayes 0.9.3

* Published on CRAN.
* Version bump in response to CRAN resubmission feedback.

# glmbayes 0.9.2

* Version bump in preparation for resubmission incorporating CRAN review feedback.

# glmbayes 0.9.1

* Wrapped OpenCL-dependent examples in `\donttest{}` for CRAN compliance.
* Reduced iteration counts in rlmb Gibbs sampler example to stay within
  CRAN example time limits on slower check machines.

# glmbayes 0.9.0

First CRAN submission. This release is a stable pre-release with a
near-complete feature set relative to earlier development builds.

## Highlights

### Bayesian Generalized Linear (glmb) and Linear (lmb) modeling functions:

  `glmb()` is a Bayesian analog for the classical `glm()` function while
  `lmb()` covers Gaussian models. Calls largely mirror those for the 
  classical functions but leverage pfamilies for prior specifications.
  Method functions largely mirror those for the classical functions. 
  Samples generated by the functions are largely iid samples 
  (no MCMC convergence dignostics are needed).

### Implemented Likelihood families/ link functions:
   
  Most of the families implemented in the `glm()` function are also implemented 
  in the `glmb()` function (the `lmb()` function covers only gaussian() families). 
  Link functions that lead to log-concave likelihood functions are generally 
  implemented.  Specifically, we have the following:
  
  **Supported likelihoods:** gaussian (identity), Poisson / quasi-Poisson
  (log), binomial / quasi-binomial (logit, probit, cloglog), Gamma (log).

### Prior Family functions:

 `pfamily` constructors are used to specify priors and play the same
  kind of role for the prior specifications as `family` constructors 
  and `link` functions play for the likelihoods. Specifically, we
  have the following:

  **Supported Priors:** Normal (all families/links), Normal–Gamma and 
  independent Normal–Gamma (gaussian families), and Gamma-on-precision 
  (gaussian and Gamma families).
  
### Prior_Setup function:
 
  The package comes with a convenient `Prior_Setup()` function that provides 
  default prior input parameters for each of the implemented models. Basic calls
  (without tailoring) mirror traditional calls to the `glmb()` and `lmb()`
  functions respectively and only require the user to provide the model formula
  and (if not the gaussian family) the family/link function. 
  
  The function can also be used to easily adjust prior specifications 
  (see documentation for details).
  
### Extensive Method functions:
  
  The package comes with extensive method functions that mirror those 
  for the classical functions.  These include dedicated `print()`,
  `summary()`, `predict()` and `simulate()` functions.

### Lower Level Modeling functions:

  The package comes with lower level modeling/simulation functions
  that advanced users can use to implement block Gibbs samplers. These
  generally come with less overhead than the `glmb()` and `lmb()` functions 
  and are called internally by the the higher level modeling functions.

### RcppParallel and OpenCL GPU Acceleration Implementations
  
  Some of the simulation functions comes with use_parallel and use_opencl options
  that speed up simulation for higher dimensional models.
  
### Extensive help files, vignettes, examples and demos

  The package also comes with extensive help files for the varios functions 
  that are complemented with a rich set of vignettes. A large number of 
  examples and demos are also availabel (see the READM.md file for a sample).

---

## Earlier development history (0.1.x series)

The notes below summarize major work during the initial development series
before the 0.9.0 pre-release.

### OpenCL and GPU acceleration

- Completed the OpenCL-based grid construction framework for large models.
- Added GPU-aware envelope sizing and improved OpenCL failure handling.
- Introduced diagnostic utilities to assess OpenCL availability and
  performance.
- Improved configure scripts to detect OpenCL and provide informative
  messages.
- Expanded OpenCL documentation and added a dedicated vignette chapter.

### Parallel CPU sampling (RcppParallel)

- Enabled parallel envelope construction and parallel iid sampling.
- Added pilot functions for large-dimension grid estimation.
- Implemented thread-safe parallel sampling for independent normal-gamma
  models.

### Core statistical improvements

- Migrated to an improved independent normal-gamma simulation algorithm.
- Added theoretical derivations for independent normal-gamma regression.
- Improved UB2 and RSS minimization routines, including scaling corrections.
- Enhanced `Prior_Setup()` to support family-specific prior construction.
- Added dedicated envelope evaluation and sizing functions.

### Package infrastructure

- Significant cleanup to remove NOTES and improve CRAN readiness.
- Improved configure and Makevars files for portability.
- Added testthat tests, including OpenCL-specific tests.
- Consolidated envelope-building functions into a cleaner structure.

### Documentation

- Major updates to README and package-level documentation.
- Added multiple new vignettes and expanded existing ones.
- Improved examples for `lmb()`, `rlmb()`, and OpenCL models.

### Bug fixes (0.1.x era)

- Corrected scaling in UB2 minimization.
- Improved error handling for missing OpenCL functionality.
- Fixed various small issues uncovered during parallelization work.
