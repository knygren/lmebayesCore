# **glmbayesCore**-only exports (by function type)

Categorized catalog of **glmbayesCore** symbols that are **not** re-exported
from **lmebayes**, organized by the kind of work each group does. Use this when
reviewing whether an export should stay public, become `@noRd`, or move with
**glmbayes** phase-out.

**Companion (overlap matrix):** [R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md)
— which symbols are shared with **glmbayes** / **lmebayes**, direct vs indirect
**lmebayes** callers, and S3 methods.

**Index:** [R_FUNCTION_INVENTORY.md](R_FUNCTION_INVENTORY.md).

---

## Scope

| Label | Meaning |
|-------|---------|
| **Core-only today** | In **glmbayesCore** `NAMESPACE`; not exported from **glmbayes** or **lmebayes** |
| **Phase-out (glmbayes)** | Still exported from **glmbayes** today; maintainer policy is **glmbayesCore**-only — call `glmbayesCore::…` or rely on `pfamily$simfun` / `Prior_Setup()` routing |
| **lmebayes direct** | **lmebayes** `R/` names `glmbayesCore::…` or `importFrom` — must stay exported in Core |
| **lmebayes indirect** | **lmebayes** reaches only via re-exported `rlmerb()` / `rglmerb()` — export optional for **lmebayes** |

Categories below follow the **glmbayes** vignette arc (simulation → envelopes →
accept–reject; Part 5 / Chapters A05–A08) and the **glmbayesCore** README split
(C++ engine, iid `simfunction`s, two-block mixed models, block-Gibbs extensions).

---

## 1. Prior calibration (internal)

Gaussian hyperparameter construction used inside default prior setup — not a
standalone user entry point (see **glmbayes** Chapter A12).

| Function | File | Status |
|----------|------|--------|
| `compute_gaussian_prior()` | `compute_gaussian_prior.R` | Phase-out (**glmbayes**); used only inside `Prior_Setup()` |

---

## 2. Simulation registry and family pipeline

Introspection and GLM family closure bundles (`f1`–`f4`, `f7`) for the envelope
path (Chapter A03 / A05).

| Function | File | Status |
|----------|------|--------|
| `simfunction()` | `simfunction.R` | Phase-out (**glmbayes**) |
| `glmbfamfunc()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |

---

## 3. Iid regression posterior samplers (`simfunction` layer)

Matrix-input draws for a single response column — the building blocks embedded
in `pfamily$simfun` and called from `rlmb()` / `rglmb()` (Chapters A05–A07).
Gaussian families use conjugate or near-conjugate paths; non-Gaussian GLMs route
to the envelope pipeline (§5–§6).

| Function | File | Status |
|----------|------|--------|
| `rNormal_reg()` | `simfunction.R` | Phase-out (**glmbayes**) |
| `rNormalGamma_reg()` | `simfunction.R` | Phase-out (**glmbayes**) |
| `rindepNormalGamma_reg()` | `simfunction.R` | Phase-out (**glmbayes**) |
| `rGamma_reg()` | `simfunction.R` | Phase-out (**glmbayes**) |
| `rGamma_Conjugate_reg()` | `simfunction.R` | Phase-out (**glmbayes**) |
| `rBeta_reg()` | `simfunction.R` | Phase-out (**glmbayes**) |

---

## 4. Standardized GLM envelope samplers and fitter hooks

Standardize the design (`glmb_Standardize_Model`), build subgradient envelopes,
and draw via `rNormalGLM_std()` / `rIndepNormalGammaReg_std()`. Weighted
working fits support mode finding inside the pipeline.

| Function | File | Status |
|----------|------|--------|
| `rNormalGLM_std()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `rIndepNormalGammaReg_std()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `glmb_Standardize_Model()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `glmb.wfit()` | `fitter_functions.R` | Phase-out (**glmbayes**) |
| `rNormal_reg.wfit()` | `fitter_functions.R` | Phase-out (**glmbayes**) |

---

## 5. Envelope construction (accept–reject geometry)

Piecewise-exponential envelope grid: size, build, evaluate, optimize, sort, and
orchestrate (Chapter A08; GPU path via OpenCL in `src/`).

| Function | File | Status |
|----------|------|--------|
| `EnvelopeSize()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `EnvelopeBuild()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `EnvelopeEval()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `EnvelopeOpt()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `EnvelopeSort()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `EnvelopeSetGrid()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `EnvelopeSetLogP()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `EnvelopeDispersionBuild()` | `simulationpipeline.R` | Phase-out (**glmbayes**) |
| `EnvelopeCentering()` | `envelopeorchestrator.R` | Phase-out (**glmbayes**) |
| `EnvelopeOrchestrator()` | `envelopeorchestrator.R` | Phase-out (**glmbayes**) |

---

## 6. Truncated distributions (C++ / OpenCL callbacks)

Low-level cdf / quantile / sampler hooks for truncated Normal, Gamma, and
inverse-Gamma draws inside ING and envelope updates (Chapters A06–A07).

| Function | File | Status |
|----------|------|--------|
| `pnorm_ct()` | `normal_ct.R` | Phase-out (**glmbayes**) |
| `rnorm_ct()` | `normal_ct.R` | Phase-out (**glmbayes**) |
| `pinvgamma_ct()` | `invgamma_ct.R` | Phase-out (**glmbayes**) |
| `qinvgamma_ct()` | `invgamma_ct.R` | Phase-out (**glmbayes**) |
| `rinvgamma_ct()` | `invgamma_ct.R` | Phase-out (**glmbayes**) |
| `rgamma_ct()` | `gamma_ct.R` | Phase-out (**glmbayes**) |

---

## 7. Row-block partition and BY-style block samplers

SAS **BY**-style row splits: normalize block indices, then run Block~1 Normal or
GLM draws conditional on hyperparameters (`inst/BLOCK_GIBBS_ERGODICITY.md`).
**lmebayes** `lmbBlock()` / `glmbBlock()` use `normalize_block()` directly;
per-block fits call **glmbayes** `lmb()` / `glmb()`, not the matrix APIs below.

| Function | File | Status | **lmebayes** |
|----------|------|--------|--------------|
| `normalize_block()` | `simfunction_block_utils.R` | Core-only today | Direct — `lmbBlock()`, `glmbBlock()`, `Prior_SetupBlock()`, `block_check_identifiability_xy()` |
| `block_rNormalReg()` | `simfunction_block.R` | Core-only today | — |
| `block_rNormalGLM()` | `simfunction_block.R` | Core-only today | — |
| `block_rNormalReg_update()` | `simfunction_block.R` | Core-only today | — |
| `block_rNormalGLM_update()` | `simfunction_block.R` | Core-only today | — |

---

## 8. Multi-response matrix samplers (`cbind` LHS)

One draw updates several response columns jointly (contrast with §7 row blocks).
Typical workflow: `multi_prior_setup()` → `multi_rlmb()` (both planned **glmbayes**
retain re-exports). **lmebayes** row blocks (`lmbBlock()`) use per-column
**glmbayes** `lmb()`, not this API.

| Function | File | Status |
|----------|------|--------|
| `multi_rlmb()` | `multi_rlmb.R` | **Planned glmbayes** retain re-export (calls `rlmb()` per column) |
| `multi_rNormalGamma_reg()` | `multi_rlmb.R` | Live via planned `multi_rlmb` API → `rNormalGamma_reg()` per column |
| `multi_rindepNormalGamma_reg()` | `multi_rlmb.R` | Live via planned `multi_rlmb` API → `rindepNormalGamma_reg()` per column |
| `multi_rNormal_reg()` | `multi_rNormal_reg.R` | Core-only today (parallel matrix API; not called by `multi_rlmb()`) |

---

## 9. Mixed-model Block~1 prep and Block~2 starts

Observation-level prior means and Gaussian ICM / GLMM mode starts for formula
drivers. Required exports for **lmebayes** (`importFrom` or qualified calls).

| Function | File | Status | **lmebayes** callers |
|----------|------|--------|----------------------|
| `build_mu_all()` | `build_mu_all.R` | Core-only today | `lmerb()`, `glmerb()` (`simulate = FALSE` → `fixef.mu`) |
| `lmerb_posterior_mean()` | `lmebayes_posterior_icm.R` | Core-only today | `lmerb()` (`simulate = FALSE`; ICM at **fixed** variance components) |
| `glmerb_posterior_mode()` | `lmebayes_posterior_icm.R` | Core-only today | `glmerb()` (`simulate = FALSE`; same) |

When `simulate = TRUE`, re-exported `rlmerb()` / `rglmerb()` perform the same
prep internally.

---

## 10. LMM / GLMM replicate-chain routers (two-block outer stage)

Replicate independent chains; each chain runs inner Gibbs sweeps between
Block~1 (random effects) and Block~2 (variance components / hyperparameters).
**lmebayes** `lmerb()` / `glmerb()` call re-exported `rlmerb()` / `rglmerb()`,
which route here — **lmebayes** never names these symbols in `R/`.

All matrix LMM engines live in **`R/rLMM_reg.R`** (shared help: **`?rLMM_reg`**).
GLMM engines live in **`R/rGLMM_reg.R`** (shared help: **`?rGLMM_reg`**).
There is no standalone **`rLMM()`** function.

| Function | File | Status | Route from **lmebayes** |
|----------|------|--------|---------------------------|
| `rLMMNormal_reg_known_vcov()` | `rLMM_reg.R` | Exported | `rlmerb()` → `.lmebayes_run_lmm_engine()` when fixed σ², all `dNormal` Block~2 |
| `rLMMNormal_reg_estimated_vcov()` | `rLMM_reg.R` | Exported | Same when fixed σ², ING Block~2 |
| `rLMMindepNormalGamma_reg_known_vcov()` | `rLMM_reg.R` | Exported | Same when dGamma `dispersion_ranef`, all `dNormal` Block~2 |
| `rLMMindepNormalGamma_reg_estimated_vcov()` | `rLMM_reg.R` | Exported | Same when dGamma σ² and ING Block~2 |
| `rLMMNormal_reg()` | `rLMM_reg.R` | Exported | Dispatcher for fixed σ² routes |
| `rLMMindepNormalGamma_reg()` | `rLMM_reg.R` | Exported | Legacy outer loop (not default **`rlmerb()`** path) |
| `rGLMM_reg()` | `rGLMM_reg.R` | Exported | `glmerb()` → `rglmerb()` when non-Gaussian `family` |
| `rGLMM_reg_known_vcov()` | `rGLMM_reg.R` | Exported | Same; all Block~2 `dNormal`; non-Gaussian always pilots; standard eigenvalue bounds |
| `rGLMM_reg_estimated_vcov()` | `rGLMM_reg.R` | Exported | Same; ING Block~2; non-Gaussian always pilots; `disp_lower` eigenvalue bounds |
| `rGLMM_sweep()` | `rGLMM_sweep.R` | Exported | Driver behind **`rGLMM_reg`** and ING LMM routes |
| `rGLMM_Re_Draw()` | `two_block_batch_gibbs.R` | Exported | Single sweep-outer re-draw helper |

See `inst/ARCHITECTURE_glmerb.md` for the full Block~1 / Block~2 call graph.

---

## 11. Two-block Gibbs engines (inner Block 1 ↔ Block 2)

Current v2 driver for cycling fixed hyperparameters and random effects with
`pfamily_list` Block~2 (used inside §10 routers and `data-raw/` calibration).

| Function | File | Status |
|----------|------|--------|
| `two_block_rNormal_reg()` | `two_block_rNormal_reg.R` | Core-only today (current R path) |

---

## 12. Two-block ergodicity and total-variation calibration

Mode weights, convergence rates, and TV bounds for inner Gibbs sweeps (Nygren
2020; `inst/BLOCK_GIBBS_ERGODICITY.md`). Used when deriving `m_convergence` and
pilot/main allocation for `rlmerb()` / `rglmerb()`.

| Function | File | Status |
|----------|------|--------|
| `two_block_rate()` | `two_block_ergodicity.R` | Core-only today |
| `two_block_rate_from_pfamily_list()` | `two_block_ergodicity.R` | Core-only today |
| `two_block_tv_bound()` | `two_block_ergodicity.R` | Core-only today |
| `two_block_l_for_tv()` | `two_block_ergodicity.R` | Core-only today |

Internal (same file): `two_block_mode_weights()` — IRLS/Fisher weights at the posterior mode for non-Gaussian rate calibration via **`rGLMM_reg`**.

---

## 13. Pilot vs main chain cost optimization

Model sampling cost vs inner-sweep count when calibrating pilot chains for ING
Block~2 and non-Gaussian GLMM paths.

| Function | File | Status |
|----------|------|--------|
| `two_block_pilot_sampling_cost()` | `two_block_pilot_cost.R` | Core-only today |
| `two_block_optimize_pilot_cost()` | `two_block_pilot_cost.R` | Core-only today |
| `two_block_d0_pilot_start()` | `two_block_pilot_cost.R` | Core-only today |
| `two_block_m_convergence_for_pilot_start()` | `two_block_pilot_cost.R` | Core-only today |

---

## 14. Two-block C++ batch path (align, Block~2 one-chain)

R and `.Call` exports for aligning random effects to hyper rows and running one
Block~2 Gibbs chain in C++ (v5 / sweep-outer batch updates).

| Function | File | Status |
|----------|------|--------|
| `two_block_align_b_to_xhyper()` | `two_block_batch_gibbs.R` | Core-only today |
| `two_block_align_b_to_xhyper_cpp()` | `two_block_batch_gibbs.R` | Core-only today |
| `two_block_block2_one_chain()` | `two_block_batch_gibbs.R` | Core-only today |
| `two_block_block2_one_chain_cpp()` | `two_block_batch_gibbs.R` | Core-only today |
| `.two_block_align_b_to_xhyper` | `two_block_batch_gibbs.R` | Core-only today (namespace alias) |
| `.two_block_block2_one_chain` | `two_block_batch_gibbs.R` | Core-only today (namespace alias) |

---

## 15. Build and compile-time diagnostics

| Function | File | Status |
|----------|------|--------|
| `glmbayesCore_has_opencl()` | `gpu_diagnostics.R` | Core-only today (distinct from **glmbayes** / **lmebayes** `has_opencl()`) |

---

## Summary counts

| Group | Symbols | Notes |
|-------|---------|-------|
| §1–§6 Phase-out from **glmbayes** | 30 | Still on **glmbayes** `NAMESPACE` until migration |
| §7–§15 Core-only today (not on **lmebayes** export surface) | 37 | Includes 4 **lmebayes** direct + 9 **lmebayes** indirect; excludes planned **glmbayes** `multi_rlmb` retain |
| **Total catalogued** | **70** | Overlap matrix: [R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md) |

---

## Review checklist

| Priority | Item |
|----------|------|
| 1 | Do not drop §9 symbols — **lmebayes** `importFrom` / direct `glmbayesCore::` depends on them. |
| 2 | §10 indirect engines may move to `@noRd` once **glmbayes** migration and `data-raw/` call sites are audited. |
| 3 | Complete **glmbayes** phase-out (§1–§6) before removing duplicate exports from **glmbayes** `NAMESPACE`. |
| 4 | Keep category boundaries aligned with vignette chapters when adding new simfunctions or envelope steps. |
