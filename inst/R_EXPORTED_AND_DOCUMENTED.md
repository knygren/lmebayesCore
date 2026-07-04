# `R/` — exported and documented functions

Symbols defined in **`R/`** that are **exported** (`NAMESPACE`) or have a **help
page** (`man/*.Rd`). Use this list when reviewing the public API, `\usage`
blocks, and README coverage.

Companion: [R_INTERNAL_HELPERS.md](R_INTERNAL_HELPERS.md) (`@noRd` and other
undocumented helpers).

Index: [R_FUNCTION_INVENTORY.md](R_FUNCTION_INVENTORY.md).

Exports below are grouped by overlap with **glmbayes** and **lmebayes**
(non-overlapping priority: glmbayes → lmebayes → lmebayes-only callers →
glmbayesCore-only).

---

## Also exported from **glmbayes** (shared API)

Present in **`NAMESPACE`** of both **glmbayesCore** and **glmbayes** (42
symbols). Signatures should stay aligned per package policy.

| Function | File | Role |
|----------|------|------|
| `Prior_Setup()` | `prior.R` | Default prior calibration from formula / design. |
| `Prior_Check()` | `prior.R` | Prior predictive checks for GLM priors. |
| `pfamily()` | `pfamily.R` | Prior-family constructors (`dNormal`, `dGamma`, …). |
| `dNormal()` | `pfamily.R` | Multivariate Normal prior component. |
| `dGamma()` | `pfamily.R` | Gamma prior component. |
| `dBeta()` | `pfamily.R` | Beta prior component. |
| `dNormal_Gamma()` | `pfamily.R` | Normal–Gamma prior component. |
| `dIndependent_Normal_Gamma()` | `pfamily.R` | Independent Normal–Gamma prior component. |
| `compute_gaussian_prior()` | `compute_gaussian_prior.R` | Gaussian prior from moments / regression fit. |
| `multi_prior_setup()` | `multi_prior_setup.R` | Multi-response prior setup from formula. |
| `simfunction()` | `simfunction.R` | Registry of simulation functions by family. |
| `glmbfamfunc()` | `simulationpipeline.R` | GLM family helpers for sampling pipeline. |
| `rglmb()` | `rglmb.R` | Matrix-level Bayesian GLM sampler. |
| `rlmb()` | `rlmb.R` | Matrix-level Bayesian LM sampler. |
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
| `diagnose_glmbayes()` | `gpu_diagnostics.R` | OpenCL / GPU diagnostic report. |

---

## Present in **lmebayes** (re-export or `importFrom`)

Exported or re-exported from **lmebayes** (`NAMESPACE` `export`, `reexports_glmbayesCore.R`), or imported via `importFrom(glmbayesCore, …)` for formula drivers — but **not** in **glmbayes** (8 symbols).

| Function | File | Role |
|----------|------|------|
| `model_setup()` | `model_setup.R` | Parse lme4-style formula → design object (`"model_setup"`). |
| `Prior_Setup_lmebayes()` | `Prior_Setup_lmebayes.R` | Calibrate Block~2 hyperpriors from reference `lmer` / `glmer`. |
| `rlmerb()` | `rlmerb.R` | Matrix-level Gaussian LMM two-block sampler. |
| `rglmerb()` | `rglmerb.R` | Matrix-level GLMM two-block sampler. |
| `pfamily_list()` | `pfamily_list.R`, `pfamily_list_lmebayes_prior_setup.R` | S3 generic; `lmebayes_prior_setup` method builds Block~2 priors. |
| `plot_sweep_history_diag()` | `plot_sweep_history_diag.R` | Cross-chain mean/SD vs inner sweep for `two_block_sweep_history`. |
| `build_mu_all()` | `build_mu_all.R` | Observation-level prior means from design + fixef (`importFrom` only). |
| `lmerb_posterior_mean()` | `lmerb_posterior_mean.R` | Gaussian LMM ICM / posterior mean for Block~2 start (`importFrom` only). |
| `glmerb_posterior_mode()` | `glmerb_posterior_mode.R` | GLMM posterior mode for Block~2 start (`importFrom` only). |

**Note:** **lmebayes** re-exports `Prior_Setup()`, `dNormal()`, `dNormal_Gamma()`,
`dIndependent_Normal_Gamma()`, and `dGamma()` from this package (listed under **glmbayes**
shared API above). S3 `print.model_setup` and `print.lmebayes_prior_setup` register in Core;
**lmebayes** `import(glmbayesCore)` dispatches them for re-exported objects.

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

**lmebayes** calls `extract_mer_variance_components()` from `summary.lmerb()` via `glmbayesCore:::`.

---

## Called from **lmebayes** / **glmbayes** without being exported there

Referenced from **lmebayes** `R/` as `glmbayesCore::…` / `glmbayesCore:::…` (or `getFromNamespace`),
exported from **glmbayesCore** only (6 symbols). **glmbayes** does not depend on **glmbayesCore**.

| Function | File | Role |
|----------|------|------|
| `rGLMM()` | `rGLMM.R` | GLMM sweep-outer engine (`rGLMM_sweep` driver). |
| `rLMMNormal_reg()` | `rLMMNormal_reg.R` | Gaussian LMM replicate-chain router. |
| `rLMMNormal_reg_estimated_vcov()` | `rLMMNormal_reg.R` | LMM with estimated residual variance. |
| `rLMMindepNormalGamma_reg()` | `rLMMNormal_reg.R` | Gaussian LMM with ING measurement dispersion. |
| `normalize_block()` | `simfunction_block_utils.R` | Row-block index normalization for block fits. |
| `multi_rlmb()` | `multi_rlmb.R` | Multi-response block LM sampler. |

---

## **glmbayesCore**-only exports

Not in **glmbayes** or **lmebayes** export surfaces and not called directly
from **lmebayes** `R/` (29 symbols).

| Function | File | Role |
|----------|------|------|
| `block_rNormalReg()` | `simfunction_block.R` | Row-block Normal regression sampler. |
| `block_rNormalGLM()` | `simfunction_block.R` | Row-block GLM envelope sampler. |
| `block_rNormalReg_update()` | `rNormalReg_reg_block_update.R` | Block~1 Normal reg draw given `mu_all`. |
| `block_rNormalGLM_update()` | `rNormalGLM_reg_block_update.R` | Block~1 GLM draw given `mu_all`. |
| `multi_rNormal_reg()` | `multi_rNormal_reg.R` | Multi-response Normal regression sampler. |
| `multi_rNormalGamma_reg()` | `multi_rlmb.R` | Multi-response Normal–Gamma sampler. |
| `multi_rindepNormalGamma_reg()` | `multi_rlmb.R` | Multi-response ING sampler. |
| `rGLMM_sweep()` | `rGLMM_sweep.R` | Sweep-outer Gibbs driver for two-block GLMM. |
| `rGLMM_Re_Draw()` | `two_block_batch_gibbs.R` | Single sweep-outer re-draw helper. |
| `rLMMNormal_reg_known_vcov()` | `rLMMNormal_reg.R` | LMM with known residual variance. |
| `two_block_rNormal_reg()` | `two_block_rNormal_reg.R` | Two-block Normal regression (legacy v1). |
| `two_block_rNormal_reg_v2()` | `two_block_rNormal_reg_v2.R` | Two-block engine with `pfamily_list` Block~2. |
| `two_block_rNormal_reg_v5()` | `two_block_rNormal_reg_v5.R` | v5 two-block driver (C++ batch path). |
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

Registered methods with help pages; not counted in the export groups above.

| Method | File | Role |
|--------|------|------|
| `formula.summary.rglmb()` | `formula.summary.rglmb.R` | Formula method for `summary.rglmb`. |
| `pfamily.default()` | `pfamily.R` | Default method for `pfamily()`. |
| `pfamily_list.lmebayes_prior_setup()` | `pfamily_list_lmebayes_prior_setup.R` | Build Block~2 `pfamily_list` from prior-setup object. |
| `print.PriorSetup()` | `prior.R` | Print prior-setup object. |
| `print.lmebayes_prior_setup()` | `Prior_Setup_lmebayes.R` | Print mixed-model prior-setup object. |
| `print.model_setup()` | `model_setup.R` | Print mixed-model design object. |
| `print.glmbfamfunc()` | `simulationpipeline.R` | Print GLM family helper object. |
| `print.pfamily()` | `pfamily.R` | Print prior-family object. |
| `print.rGamma_reg()` | `simfunction.R` | Print Gamma-regression sample object. |
| `print.rglmb()` / `print.rlmb()` | `rglmb.R`, `rlmb.R` | Print GLM / LM sample objects. |
| `print.simfunction()` | `simfunction.R` | Print simulation-function registry. |
| `print.summary.rGamma_reg()` | `summary.rgamma_reg.R` | Print Gamma-regression summary. |
| `print.summary.rglmb()` | `summary.rglmb.R` | Print GLM sample summary. |
| `print.summary.mrglmb()` | `summary.mrglmb.R` | Print multi-response GLM summary. |
| `print.two_block_sweep_history()` | `two_block_sweep_history.R` | Print sweep-history diagnostics. |
| `print.two_block_rate()` | `two_block_rate.R` | Print rate object with TV tolerances. |
| `print.two_block_mode_weights()` | `two_block_mode_weights.R` | Print mode-weight table. |
| `residuals.rglmb()` (+ fitted, working) | `residuals.rglmb.R` | Residuals methods for `rglmb`. |
| `residuals.rlmb()` | `residuals.rglmb.R` | Residuals for `rlmb`. |
| `simfunction.default()` | `simfunction.R` | Default simfunction method. |
| `summary.rGamma_reg()` | `summary.rgamma_reg.R` | Summary for Gamma-regression samples. |
| `summary.rglmb()` / `summary.rlmb()` | `summary.rglmb.R`, `rlmb.R` | Posterior summaries for samples. |
| `summary.mrglmb()` | `summary.mrglmb.R` | Multi-response GLM summary. |

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
| 2 | Treat **lmebayes** call surface (`rGLMM`, `normalize_block`, …) as a semver-sensitive API. |
| 3 | Run `devtools::document()` after any `@export` or `\usage` change. |
