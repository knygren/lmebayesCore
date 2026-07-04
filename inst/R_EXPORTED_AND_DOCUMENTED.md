# `R/` — exported and documented functions

Symbols defined in **`R/`** that are **exported** (`NAMESPACE`) or have a **help
page** (`man/*.Rd`). Use this list when reviewing the public API, `\usage`
blocks, and README coverage.

Companion: [R_INTERNAL_HELPERS.md](R_INTERNAL_HELPERS.md) (`@noRd` and other
undocumented helpers).

**Core-only catalog (by function type):**
[R_CORE_ONLY_EXPORTS.md](R_CORE_ONLY_EXPORTS.md).

**Reachability / dead-export analysis:**
[R_EXPORT_REACHABILITY.md](R_EXPORT_REACHABILITY.md).

Index: [R_FUNCTION_INVENTORY.md](R_FUNCTION_INVENTORY.md).

Exports below are grouped by overlap with **glmbayes** and **lmebayes**
(non-overlapping priority: glmbayes → lmebayes → lmebayes-only callers →
glmbayesCore-only).

---

## Also exported from **glmbayes** (shared API)

Present in **`NAMESPACE`** of both **glmbayesCore** and **glmbayes** (42
symbols). Signatures should stay aligned per package policy while both packages
export a symbol.

**Migration policy:** **glmbayes** should keep re-exporting the *user-facing*
prior and sampler API (`Prior_Setup`, `pfamily` constructors, `rglmb` / `rlmb`,
`multi_prior_setup`, `multi_rlmb`, …). Internal calibration helpers (e.g.
`compute_gaussian_prior`, used only inside `Prior_Setup()`), plus lower-level
simulation, envelope, and C++ callback exports, become **glmbayesCore**-only
over time — call them as `glmbayesCore::…` (or rely on `Prior_Setup()` /
`pfamily$simfun` routing) rather than expecting them on `library(glmbayes)`.

S3 methods registered in **glmbayes** today are listed separately from Core-only
methods at the end of this document (see **S3 methods**).

### 1) Retain as **glmbayes** re-exports (core API)

End-user and formula-level workflow; stay exported from **glmbayes** (13 symbols;
`multi_rlmb()` planned — not yet in **glmbayes** `NAMESPACE`).

| Function | File | Role |
|----------|------|------|
| `Prior_Setup()` | `prior.R` | Default prior calibration from formula / design. |
| `Prior_Check()` | `prior.R` | Prior predictive checks for GLM priors. |
| `pfamily()` | `pfamily.R` | Prior-family generic (`dNormal`, `dGamma`, …). |
| `dNormal()` | `pfamily.R` | Multivariate Normal prior component. |
| `dGamma()` | `pfamily.R` | Gamma prior component. |
| `dBeta()` | `pfamily.R` | Beta prior component. |
| `dNormal_Gamma()` | `pfamily.R` | Normal–Gamma prior component. |
| `dIndependent_Normal_Gamma()` | `pfamily.R` | Independent Normal–Gamma prior component. |
| `multi_prior_setup()` | `multi_prior_setup.R` | Multi-response Gaussian prior setup (`cbind` LHS). |
| `multi_rlmb()` | `multi_rlmb.R` | Multi-response LM draws (`cbind` LHS; `rlmb()` per column). Pairs with `multi_prior_setup()`. |
| `rglmb()` | `rglmb.R` | Matrix-level Bayesian GLM sampler (`glmb()` backend). |
| `rlmb()` | `rlmb.R` | Matrix-level Bayesian LM sampler (`lmb()` backend). |
| `diagnose_glmbayes()` | `gpu_diagnostics.R` | OpenCL / GPU diagnostic report. |

**Co-exports (same help topic, planned with `multi_rlmb()`):**
`multi_rNormalGamma_reg()`, `multi_rindepNormalGamma_reg()`.

### 2) Phase out of **glmbayes** (low-level; **glmbayesCore**-only)

Internal calibration, simulation registry, envelope pipeline, fitter hooks, and
truncated-distribution C++ callbacks (30 symbols). Still exported from
**glmbayesCore**; remove from **glmbayes** `NAMESPACE` when downstream call
sites use `glmbayesCore::` or internal routing only.

| Function | File | Role |
|----------|------|------|
| `compute_gaussian_prior()` | `compute_gaussian_prior.R` | Gaussian calibration inside `Prior_Setup()` only (not a user entry point). |
| `simfunction()` | `simfunction.R` | Registry of simulation functions by family. |
| `glmbfamfunc()` | `simulationpipeline.R` | GLM family helpers for sampling pipeline. |
| `rNormal_reg()` | `simfunction.R` | Normal regression posterior sampler. |
| `rNormal_reg.wfit()` | `fitter_functions.R` | Weighted Normal regression fit helper. |
| `rNormalGamma_reg()` | `simfunction.R` | Normal–Gamma regression sampler. |
| `rindepNormalGamma_reg()` | `simfunction.R` | Independent Normal–Gamma regression sampler. |
| `rGamma_reg()` | `simfunction.R` | Gamma regression sampler. |
| `rGamma_Conjugate_reg()` | `simfunction.R` | Conjugate Gamma regression sampler. |
| `rBeta_reg()` | `simfunction.R` | Beta regression sampler. |
| `rNormalGLM_std()` | `simulationpipeline.R` | Standardized Normal–GLM envelope sampler. |
| `rIndepNormalGammaReg_std()` | `simulationpipeline.R` | Standardized ING GLM envelope sampler. |
| `glmb.wfit()` | `fitter_functions.R` | Weighted GLM working fit for sampling. |
| `glmb_Standardize_Model()` | `simulationpipeline.R` | Standardize design for envelope sampling. |
| `EnvelopeSize()` | `simulationpipeline.R` | Envelope grid size from subgradient geometry. |
| `EnvelopeBuild()` | `simulationpipeline.R` | Build accept–reject envelope. |
| `EnvelopeEval()` | `simulationpipeline.R` | Evaluate envelope log-density. |
| `EnvelopeOpt()` | `simulationpipeline.R` | Optimize envelope grid / core count. |
| `EnvelopeSort()` | `simulationpipeline.R` | Sort envelope grid for sampling. |
| `EnvelopeSetGrid()` | `simulationpipeline.R` | Assign envelope grid indices. |
| `EnvelopeSetLogP()` | `simulationpipeline.R` | Set envelope log-probability table. |
| `EnvelopeDispersionBuild()` | `simulationpipeline.R` | Dispersion envelope for ING models. |
| `EnvelopeCentering()` | `envelopeorchestrator.R` | Center envelope at posterior mode. |
| `EnvelopeOrchestrator()` | `envelopeorchestrator.R` | Full envelope build + dispersion orchestration. |
| `pnorm_ct()` / `rnorm_ct()` | `normal_ct.R` | Truncated Normal cdf / sampler (C++/OpenCL). |
| `pinvgamma_ct()` / `qinvgamma_ct()` / `rinvgamma_ct()` | `invgamma_ct.R` | Inverse-Gamma cdf / quantile / sampler. |
| `rgamma_ct()` | `gamma_ct.R` | Truncated Gamma sampler. |

---

## Present in **lmebayes** (from **glmbayesCore**)

Mixed-model symbols defined in **glmbayesCore** and used by **lmebayes**, but
**not** in **glmbayes** (6 re-exported; 4 additional direct Core calls — see below).

**Note:** **lmebayes** also re-exports `Prior_Setup()`, `dNormal()`,
`dNormal_Gamma()`, `dIndependent_Normal_Gamma()`, and `dGamma()` (listed under
**glmbayes** shared API above). S3 `print.model_setup`,
`print.lmebayes_prior_setup`, and `pfamily_list.lmebayes_prior_setup` register
in Core; **lmebayes** `import(glmbayesCore)` dispatches them for re-exported
objects.

### Re-exports (6 symbols)

Exported from **lmebayes** via `R/reexports_glmbayesCore.R` (`@export` /
`importFrom` + re-export).

| Function | File | Role |
|----------|------|------|
| `model_setup()` | `model_setup.R` | Parse lme4-style formula → design object (`"model_setup"`). |
| `Prior_Setup_lmebayes()` | `Prior_Setup_lmebayes.R` | Calibrate Block~2 hyperpriors from reference `lmer` / `glmer`. |
| `rlmerb()` | `rlmerb.R` | Matrix-level Gaussian LMM two-block sampler. |
| `rglmerb()` | `rglmerb.R` | Matrix-level GLMM two-block sampler. |
| `pfamily_list()` | `pfamily_list.R`, `pfamily_list_lmebayes_prior_setup.R` | S3 generic; `lmebayes_prior_setup` method builds Block~2 priors. |
| `plot_sweep_history_diag()` | `plot_sweep_history_diag.R` | Cross-chain mean/SD vs inner sweep for `two_block_sweep_history`. |

**lmebayes** callers: `model_setup()` and `Prior_Setup_lmebayes()` from user
workflows and `Prior_SetupBlock()`; `pfamily_list()` from `lmerb()` / `glmerb()`
(via `.lmebayes_priors_from_pfamily_list()`); `rlmerb()` / `rglmerb()` from
`lmerb()` / `glmerb()` when `simulate = TRUE`; `plot_sweep_history_diag()` from
demos and user diagnostics on `fit$sweep_history$main`.

### Direct calls from **lmebayes** `R/` (must stay exported in Core)

**lmebayes** references these as `glmbayesCore::…` (or `importFrom`). They must
remain in **glmbayesCore** `NAMESPACE` for **lmebayes** to load.

| Function | File | **lmebayes** callers | **lmebayes** surface |
|----------|------|----------------------|----------------------|
| `build_mu_all()` | `build_mu_all.R` | `lmerb()`, `glmerb()` | `importFrom` only |
| `lmerb_posterior_mean()` | `lmerb_posterior_mean.R` | `lmerb()` | `importFrom` only |
| `glmerb_posterior_mode()` | `glmerb_posterior_mode.R` | `glmerb()` | `importFrom` only |
| `normalize_block()` | `simfunction_block_utils.R` | `lmbBlock()`, `glmbBlock()`, `Prior_SetupBlock()` (via `.blmb_formula_block_meta()`); `block_check_identifiability_xy()` | direct `glmbayesCore::` (not in **lmebayes** `NAMESPACE`) |

Used when `simulate = FALSE` for the three posterior/mean helpers; when
`simulate = TRUE`, re-exported `rlmerb()` / `rglmerb()` run the same prep
internally.

**Direct via `glmbayesCore:::`** (internal in Core — export not required):
`extract_mer_variance_components()` ← `summary.lmerb()`; `.lmebayes_*` helpers
← `lmerb()`, `glmerb()`. Listed under **Documented but not exported** below.

Indirect engines (`rGLMM`, `rLMMNormal_reg`, …) are listed under
**glmbayesCore-only exports** below — **lmebayes** reaches them only via
re-exported `rlmerb()` / `rglmerb()` (export optional for **lmebayes**).

---

## Documented but not exported (`@keywords internal` + `man/`)

Callable with `glmbayesCore:::`; have help pages but are not in `NAMESPACE`.

| Function | File | Role |
|----------|------|------|
| `is_single_factor_model()` | `lme4_design_utilities.R` | Exactly one grouping factor in formula. |
| `is_fixed_effects_only()` | `lme4_design_utilities.R` | No random-effects terms. |
| `get_lme4_components()` | `lme4_design_utilities.R` | lme4 parse → design matrices. |
| `show_lme4_Z_random()` | `lme4_design_utilities.R` | Debug random-effects design. |
| `classify_lme4_fixed_columns()` | `lme4_design_utilities.R` | Population vs group-level fixed columns. |
| `classify_crosslevel_re_moderation()` | `lme4_design_utilities.R` | Cross-level RE structure. |
| `extract_re_hyper_matrices()` | `lme4_design_utilities.R` | Block~2 group designs. |
| `extract_re_Z_obs()` | `lme4_design_utilities.R` | Obs-level Z for one group. |
| `extract_lme4_submatrices()` | `lme4_design_utilities.R` | Subset parsed lme4 parts. |
| `extract_lme4_fixed_group_matrix()` | `lme4_design_utilities.R` | Group-level fixed matrix. |
| `extract_lmer_variance_components()` | `lme4_design_utilities.R` | Variance components from **lmer** fit. |
| `extract_mer_variance_components()` | `lme4_design_utilities.R` | Variance components from **glmer** fit. |
| `lmerb_default_vcov_formula()` | `lme4_design_utilities.R` | Default vcov formula for prior scaling. |

**lmebayes** calls `extract_mer_variance_components()` from `summary.lmerb()` via
`glmbayesCore:::` (direct call; correctly **not** exported).

---

## **glmbayesCore**-only exports

Not in **glmbayes** or **lmebayes** export surfaces (29 symbols). None are
direct `glmbayesCore::…` references from **lmebayes** `R/`.

**By function type:** [R_CORE_ONLY_EXPORTS.md](R_CORE_ONLY_EXPORTS.md) (full
66-symbol catalog including **glmbayes** phase-out groups §1–§6).

### Indirect from **lmebayes** via `rlmerb()` / `rglmerb()` (4 symbols)

**lmebayes** never names these in `R/`; **glmbayesCore** routes to them inside
`rlmerb()`, `rglmerb()`, or `.lmebayes_run_lmm_engine()`. They do **not** need
to stay `@export` for **lmebayes** to work (could become `@noRd` internals).
They may remain exported for power users, `data-raw/`, or **glmbayes** migration.

| Function | File | Route from **lmebayes** | Role |
|----------|------|-------------------------|------|
| `rGLMM()` | `rGLMM.R` | `glmerb()` → `rglmerb()` when `simulate = TRUE`, non-Gaussian `family` | GLMM sweep-outer engine. |
| `rLMMNormal_reg()` | `rLMMNormal_reg.R` | `lmerb()` / `glmerb()` → `rlmerb()` / `rglmerb()` → `.lmebayes_run_lmm_engine()` when Block~2 is all `dNormal` and `dispersion_ranef` is a fixed scalar | Gaussian LMM replicate-chain router. |
| `rLMMNormal_reg_estimated_vcov()` | `rLMMNormal_reg.R` | Same chain when `prior$any_non_normal` (e.g. ING Block~2) | LMM with estimated residual variance. |
| `rLMMindepNormalGamma_reg()` | `rLMMNormal_reg.R` | Same chain when `dispersion_ranef` is `dGamma()` (`disp_info$mode == "gamma"`) | LMM with ING measurement dispersion. |

### Other **glmbayesCore**-only exports (25 symbols)

Not called from **lmebayes** `R/` (directly or via formula drivers).

| Function | File | Role |
|----------|------|------|
| `block_rNormalReg()` | `simfunction_block.R` | Row-block Normal regression sampler. |
| `block_rNormalGLM()` | `simfunction_block.R` | Row-block GLM envelope sampler. |
| `block_rNormalReg_update()` | `rNormalReg_reg_block_update.R` | Block~1 Normal reg draw given `mu_all`. |
| `block_rNormalGLM_update()` | `rNormalGLM_reg_block_update.R` | Block~1 GLM draw given `mu_all`. |
| `multi_rNormal_reg()` | `multi_rNormal_reg.R` | Multi-response Normal regression sampler (matrix API; not called by `multi_rlmb()`). |
| `rGLMM_sweep()` | `rGLMM_sweep.R` | Sweep-outer Gibbs driver for two-block GLMM. |
| `rGLMM_Re_Draw()` | `two_block_batch_gibbs.R` | Single sweep-outer re-draw helper. |
| `rLMMNormal_reg_known_vcov()` | `rLMMNormal_reg.R` | LMM with known residual variance. |
| `two_block_rNormal_reg_v2()` | `two_block_rNormal_reg_v2.R` | Two-block engine with `pfamily_list` Block~2. |
| `two_block_rate()` | `two_block_rate.R` | Block~2 convergence rate from mode weights. |
| `two_block_rate_v2()` | `two_block_rNormal_reg_v2.R` | Rate helper aligned with v2 sampler. |
| `two_block_mode_weights()` | `two_block_mode_weights.R` | Mode weights for rate / TV calibration. |
| `two_block_tv_bound()` | `two_block_tv_bound.R` | Total-variation bound vs inner sweeps. |
| `two_block_l_for_tv()` | `two_block_tv_bound.R` | Invert TV bound for target tolerance. |
| `two_block_pilot_sampling_cost()` | `two_block_pilot_cost.R` | Pilot chain sampling cost model. |
| `two_block_optimize_pilot_cost()` | `two_block_pilot_cost.R` | Optimize pilot vs main chain allocation. |
| `two_block_d0_pilot_start()` | `two_block_pilot_cost.R` | Default pilot start for cost optimization. |
| `two_block_m_convergence_for_pilot_start()` | `two_block_pilot_cost.R` | Inner sweeps for pilot-start calibration. |
| `two_block_align_b_to_xhyper()` | `two_block_batch_gibbs.R` | Align ranef vector to hyper design rows. |
| `two_block_align_b_to_xhyper_cpp()` | `two_block_batch_gibbs.R` | C++ export for align helper. |
| `two_block_block2_one_chain()` | `two_block_batch_gibbs.R` | One Block~2 Gibbs chain (fixef + tau2). |
| `two_block_block2_one_chain_cpp()` | `two_block_batch_gibbs.R` | C++ export for Block~2 one-chain draw. |
| `.two_block_align_b_to_xhyper` | `two_block_batch_gibbs.R` | Namespace alias of `two_block_align_b_to_xhyper`. |
| `.two_block_block2_one_chain` | `two_block_batch_gibbs.R` | Namespace alias of `two_block_block2_one_chain`. |
| `glmbayesCore_has_opencl()` | `gpu_diagnostics.R` | Compile-time OpenCL flag for this build. |

---

## S3 methods (`NAMESPACE` → `S3method`)

Registered in **glmbayesCore**; not counted in the export groups above. When both
**glmbayes** and **glmbayesCore** attach, **R** warns that **glmbayes** overwrites
duplicate registrations — during migration, **glmbayes** should `import(glmbayesCore)`
and stop re-registering methods that live in Core.

**Verification source:** **glmbayes** `NAMESPACE` (installed R 4.6 library) cross-checked
line-by-line against **glmbayesCore** `NAMESPACE` (**44** vs **26** registrations).
Older **glmbayes** source trees may omit some lines (e.g. `summary.mlmb`); treat the
installed namespace as authoritative.

### Complete S3 overlap (**glmbayes** ∩ **glmbayesCore**) — 15 methods

Exhaustive intersection — no other shared `S3method` registrations exist. Methods
like `extractAIC.rglmb()`, `deviance.rglmb()`, and all `*.glmb` / `*.lmb` /
`*.mlmb` dispatch are **glmbayes-only** (see below).

| Method | Core file | Policy |
|--------|-----------|--------|
| `pfamily.default()` | `pfamily.R` | **Retain** |
| `print.PriorSetup()` | `prior.R` | **Retain** |
| `print.pfamily()` | `pfamily.R` | **Retain** |
| `print.rglmb()` | `rglmb.R` | **Retain** |
| `print.rlmb()` | `rlmb.R` | **Retain** |
| `residuals.rglmb()` | `residuals.rglmb.R` | **Retain** (also `rlmb` via class) |
| `summary.rglmb()` | `summary.rglmb.R` | **Retain** (also `rlmb` via class) |
| `print.summary.rglmb()` | `summary.rglmb.R` | **Retain** |
| `formula.summary.rglmb()` | `formula.summary.rglmb.R` | **Retain** |
| `print.glmbfamfunc()` | `simulationpipeline.R` | **Phase-out** (§2) |
| `print.simfunction()` | `simfunction.R` | **Phase-out** (§2) |
| `simfunction.default()` | `simfunction.R` | **Phase-out** (§2) |
| `print.rGamma_reg()` | `simfunction.R` | **Phase-out** (§2) |
| `summary.rGamma_reg()` | `summary.rgamma_reg.R` | **Phase-out** (§2) |
| `print.summary.rGamma_reg()` | `summary.rgamma_reg.R` | **Phase-out** (§2) |

**Attach warning:** when **lmebayes** loads both packages, **R** may list a subset
of these 15 as overwritten by **glmbayes**; the table above is the full set.

### **glmbayes** S3 not in **glmbayesCore** — 29 methods

Formula-layer and **glmbayes**-specific dispatch; stay on **glmbayes** (do not
duplicate in Core).

| Method | Typical object | Notes |
|--------|----------------|-------|
| `anova.glmb()` | `glmb` | |
| `case.names.glmb()` | `glmb` | |
| `confint.glmb()` | `glmb` | |
| `cooks.distance.glmb()` | `glmb` | |
| `dfbetas.glmb()` | `glmb` | |
| `dummy.coef.glmb()` | `glmb` | |
| `extractAIC.glmb()` | `glmb` | |
| `influence.glmb()` | `glmb` | |
| `logLik.glmb()` | `glmb` | |
| `plot.glmb()` | `glmb` | |
| `predict.glmb()` | `glmb` | |
| `residuals.glmb()` | `glmb` | |
| `rstandard.glmb()` | `glmb` | |
| `rstudent.glmb()` | `glmb` | |
| `simulate.glmb()` | `glmb` | |
| `variable.names.glmb()` | `glmb` | |
| `vcov.glmb()` | `glmb` | |
| `print.glmb()` | `glmb` | |
| `print.dummy_coef.glmb()` | `dummy_coef.glmb` | |
| `print.lmb()` | `lmb` | |
| `print.mlmb()` | `mlmb` | |
| `print.directional_tail()` | `directional_tail` | |
| `residuals.lmb()` | `lmb` | |
| `deviance.rglmb()` | `rglmb` | **Not** in Core (Core has no `deviance.rglmb`) |
| `extractAIC.rglmb()` | `rglmb` | **Not** in Core |
| `summary.glmb()` | `glmb` | |
| `print.summary.glmb()` | `summary.glmb` | |
| `summary.mlmb()` | `mlmb` | Installed **glmbayes** only (may be absent in older source) |
| `print.summary.mlmb()` | `summary.mlmb` | Same |

### Summary methods — complete **glmbayes** inventory

Every `summary`, `print.summary`, and `formula.summary` registration in
**glmbayes** `NAMESPACE`. None omitted.

| Method | **glmbayesCore**? | Policy |
|--------|-------------------|--------|
| `summary.glmb()` | No | **glmbayes-only** — formula `glmb()` / `lmb()` fit objects; stays on **glmbayes** (not in Core). |
| `print.summary.glmb()` | No | **glmbayes-only** — print method for `summary.glmb`. |
| `summary.mlmb()` | No | **glmbayes-only** — formula `mlmb()` multi-response fits; stays on **glmbayes**. |
| `print.summary.mlmb()` | No | **glmbayes-only** — print method for `summary.mlmb`. |
| `summary.rglmb()` | Yes (`summary.rglmb.R`) | **Overlap / retain** — matrix sampler `rglmb()`; also serves `rlmb` via class `c("rlmb","rglmb",…)`. |
| `print.summary.rglmb()` | Yes (`summary.rglmb.R`) | **Overlap / retain** |
| `formula.summary.rglmb()` | Yes (`formula.summary.rglmb.R`) | **Overlap / retain** |
| `summary.rGamma_reg()` | Yes (`summary.rgamma_reg.R`) | **Overlap / phase-out** — drop from **glmbayes** when `rGamma_reg()` moves Core-only (§2). |
| `print.summary.rGamma_reg()` | Yes (`summary.rgamma_reg.R`) | **Overlap / phase-out** |

**Summary methods in glmbayesCore but not glmbayes** (do not add to **glmbayes**
`NAMESPACE` unless policy changes):

| Method | File | Policy |
|--------|------|--------|
| `summary.rlmb()` | `summary.rglmb.R` | **Core alias** — delegates to `summary.rglmb()`; **glmbayes** uses `summary.rglmb` for `rlmb` objects. |
| `summary.mrglmb()` | `summary.mrglmb.R` | **Retain gap** — pairs with planned `multi_rlmb()` re-export. |
| `print.summary.mrglmb()` | `summary.mrglmb.R` | **Retain gap** |

### Other S3 — overlap with **glmbayes** (non-summary) — 10 methods

Duplicate registration in both packages (excluding summary rows above).

| Method | File | Role |
|--------|------|------|
| `pfamily.default()` | `pfamily.R` | Default method for retain generic `pfamily()`. |
| `print.PriorSetup()` | `prior.R` | Print prior-setup object. |
| `print.glmbfamfunc()` | `simulationpipeline.R` | Print GLM family helper object. |
| `print.pfamily()` | `pfamily.R` | Print prior-family object. |
| `print.rGamma_reg()` | `simfunction.R` | Print Gamma-regression sample object. |
| `print.rglmb()` | `rglmb.R` | Print GLM sample object. |
| `print.rlmb()` | `rlmb.R` | Print LM sample object. |
| `print.simfunction()` | `simfunction.R` | Print simulation-function registry. |
| `residuals.rglmb()` (+ fitted, working) | `residuals.rglmb.R` | Residuals for `rglmb` and (via class) `rlmb` objects. |
| `simfunction.default()` | `simfunction.R` | Default method for `simfunction()`. |

**Phase-out note (non-summary overlap):** `print.glmbfamfunc`, `print.simfunction`,
`simfunction.default`, and `print.rGamma_reg` dispatch on phase-out exports (§2).

### Not in **glmbayes** (Core-only dispatch) — 9 methods

Registered only in **glmbayesCore**; do **not** add to **glmbayes** `NAMESPACE`
(excluding summary rows already listed above).

| Method | File | Subgroup | Role |
|--------|------|----------|------|
| `residuals.rlmb()` | `residuals.rglmb.R` | Core alias | Same draw logic as `residuals.rglmb()`; **glmbayes** uses `residuals.rglmb` for `rlmb`. |
| `pfamily_list.lmebayes_prior_setup()` | `pfamily_list_lmebayes_prior_setup.R` | lmebayes | Build Block~2 `pfamily_list` from prior-setup object. |
| `print.lmebayes_prior_setup()` | `Prior_Setup_lmebayes.R` | lmebayes | Print mixed-model prior-setup object. |
| `print.model_setup()` | `model_setup.R` | lmebayes | Print mixed-model design object. |
| `print.two_block_sweep_history()` | `two_block_sweep_history.R` | lmebayes | Print sweep-history diagnostics (`rglmerb()` / `rlmerb()`). |
| `print.two_block_rate()` | `two_block_rate.R` | Core-only | Print rate object with TV tolerances. |
| `print.two_block_mode_weights()` | `two_block_mode_weights.R` | Core-only | Print mode-weight table. |
| `residuals.summary.rglmb()` | `residuals.rglmb.R` | Core-only | Residuals for `summary.rglmb` object. |

**S3 counts:** **glmbayes** 44; **glmbayesCore** 26; **overlap** 15; **glmbayes-only** 29;
**Core-only** 11 (26 − 15).

---

## Documentation topics (no function body in `R/`)

| Topic / file | Contents |
|--------------|----------|
| `glmbayesCore-package.R` | Package meta, imports (`"_PACKAGE"`). |
| `gpu_diagnostics.R` | Links to OpenCL diagnostics topics. |
| `data-*.R` | Lazy data docs (`Boston_centered`, `BikeSharing`, `carinsca`, etc.). |

---

## Review checklist (exports / docs)

| Priority | Item |
|----------|------|
| 1 | Keep **glmbayes**-shared exports signature-aligned when touching `R/`. |
| 2 | Trim **glmbayes** re-exports in group **2** (simfunctions, envelopes, `*_ct`) once call sites use `glmbayesCore::`. |
| 3 | Treat **lmebayes** direct call surface (`build_mu_all`, `normalize_block`, …) as semver-sensitive; indirect engines (`rGLMM`, `rLMMNormal_reg`, …) may be internalized. |
| 4 | Run `devtools::document()` after any `@export` or `\usage` change. |
| 5 | When **glmbayes** imports **glmbayesCore**, drop duplicate S3 registrations (see **S3 methods** — overlap tables; keep **glmbayes-only** summary methods for `glmb` / `mlmb`). |
