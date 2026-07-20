# lmebayesCore

![GitHub release (latest by date)](https://img.shields.io/github/v/release/knygren/lmebayesCore?label=version)
![License: GPL-2](https://img.shields.io/badge/license-GPL--2-blue.svg)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/knygren/lmebayesCore/R-CMD-check.yaml?label=R%20CMD%20Check)

**lmebayesCore** is the compiled two-block Gibbs sampling engine that powers [lmebayes](https://github.com/knygren/lmebayes). It started as a history-preserving fork of [glmbayesCore](https://github.com/knygren/glmbayesCore) (2026-07-15) and has since been trimmed down to the mixed-model layer only: formula-to-design-matrix setup (`model_setup()`), Block~2 hyperprior calibration (`Prior_Setup_lmebayes()`), and the matrix-level two-block samplers (`rlmerb()`/`rglmerb()`, the `rLMM_reg`/`rGLMM_reg` route families, and row-block engines). The iid GLM/LM prior, envelope, and sampling machinery that `lmebayesCore` originally forked now lives exclusively in `glmbayesCore` (an `Imports:` dependency); `lmebayesCore` calls it via `glmbayesCore::…` rather than shipping its own copy. End users should install `lmebayes` rather than this package directly.

The relationship to the broader ecosystem parallels how `StanHeaders` / `rstan` serve as the compiled backbone for `rstanarm`: `lmebayesCore` is the infrastructure layer for mixed models; `lmebayes` is the user-facing package built on top of it. (`glmbayesCore` plays the same role for the iid `glmbayes` package, and `lmebayesCore` now depends on it directly for that iid layer.)

---

## Package Ecosystem

```
                ┌─────────────────────────────────────────┐
                │           End-user packages             │
                │        glmbayes   ·   lmebayes           │
                └──────────┬──────────────────┬────────────┘
                           │ Imports          │ Imports
                ┌──────────▼──────────┐       │
                │     glmbayesCore     │◄──────┤ Imports
                │  iid C++ samplers    │       │
                │  pfamily · priors    │       │
                │  simfunctions        │       │
                │  rglmb / rlmb        │       │
                └──────────┬───────────┘       │
                           │ Imports    ┌───────▼────────────┐
                           └────────────►     lmebayesCore     │
                                        │  two-block Gibbs     │
                                        │  mixed-model engine  │
                                        │  model_setup /       │
                                        │  rlmerb / rglmerb     │
                                        └──────────┬────────────┘
                                                   │ Imports
                                        ┌──────────▼──────────────┐
                                        │ opencltools · nmathopencl │
                                        │ Rcpp · RcppArmadillo ·    │
                                        │ RcppParallel              │
                                        └───────────────────────────┘
```

**glmbayes** provides the formula interface (`glmb()`, `lmb()`) for iid GLM/LM fitting, MCMC diagnostics, and the S3 methods that mirror base-R's `lm()` / `glm()` -- built directly on **glmbayesCore**.

**lmebayes** extends this to linear and generalized linear **mixed-effects** models (`lmerb()`, `glmerb()`), built on **lmebayesCore** for mixed-model setup/sampling and re-exporting priors from **glmbayesCore** for the iid layer inside each Gibbs block.

**Important caveat on the diagram above:** at the *R* level, `lmebayesCore` cleanly delegates all iid functionality to `glmbayesCore`. At the *C++* level, the compiled engine has **not yet** been pruned (Stage 3 of the ongoing deduplication effort) -- `lmebayesCore`'s `src/` still links its own copies of the iid envelope/sampler object code alongside the mixed-model-only `.cpp` files, because the two-block Gibbs C++ paths call directly into that iid code within the same DLL. See [Future plans](#future-plans) below.

---

## What Is Inside lmebayesCore

### C++ sampling engine (`src/`)

The core is organized under the `glmbayes::` namespace (kept as an internal implementation detail; not renamed during the fork), partitioned into sub-namespaces:

| Sub-namespace | Key files | Role |
|---|---|---|
| `glmbayes::fam` | `famfuncs.h`, `famfuncs_*.cpp` | Negative log-posterior (`f2`) and gradient (`f3`) for gaussian, poisson, binomial, Gamma |
| `glmbayes::env` | `EnvelopeBuild*.cpp`, `EnvelopeEval.cpp`, `EnvelopeSort.cpp`, `EnvelopeSize.cpp`, `Set_Grid.cpp`, `Set_LogP.cpp` | Piecewise-exponential envelope construction (Nygren & Nygren, 2006) |
| `glmbayes::sim` | `rNormalGLM.cpp`, `rIndepNormalGammaReg.cpp`, `rNormalGammaReg.cpp`, `rNormalReg.cpp`, `rGammaGamma.cpp`, `rGammaGaussian.cpp` | iid posterior samplers |
| `glmbayes::rng` | `rng_utils.cpp` | Thread-safe RNG wrappers for parallel sampling |
| `glmbayes::progress` | `progress_utils.cpp` | Optional progress bar support |

**Mixed-model-only `.cpp` files** (not present in `glmbayesCore`): `twoBlockGibbs.cpp`, `twoBlockGibbsStaged.cpp`, `two_block_block1.cpp`, `block_rIndepNormalGammaReg.cpp`, `block_utils.cpp`, `rNormalGLMBlocks.cpp`, `rNormalRegBlocks.cpp`. These implement the row-block and two-block Gibbs drivers and are the reason the iid C++ engine above is still linked into this package (Stage 3 pruning is pending -- see [Future plans](#future-plans)).

Export wrappers in `export_wrappers.cpp` and `kernel_wrappers.cpp` expose selected entry points to R via Rcpp. `src/package_ns.h` defines `GLMBAYES_R_NS` (`"glmbayesCore"`), used by `src/R_interface.h` and `src/twoBlockGibbs.cpp` to resolve R-level callbacks (`EnvelopeOpt`, `EnvelopeSort`, `glmbfamfunc`, `rNormal_reg.wfit`, `rgamma_ct`, Block~2's `rglmb`) from **glmbayesCore**'s namespace rather than a local copy.

### OpenCL kernels (`inst/cl/`)

For systems with an OpenCL-capable device, envelope construction can be offloaded to the GPU. The `inst/cl/` tree contains:

- **`src/f2_f3_*.cl`** -- GPU ports of the `f2` (negative log-posterior) and `f3` (gradient) functions for each family/link combination (gaussian, poisson, binomial logit/probit/cloglog, Gamma log).
- **`nmath/`** -- OpenCL port of the R Mathlib probability functions (`dnorm`, `dgamma`, `dbinom`, `dpois`, `lgamma`, etc.) needed by the kernels.
- **`libR_shims/`, `R_ext_*/`, `R_shims/`** -- Shim headers that make the nmath kernels compile cleanly under OpenCL's C99-based dialect.
- **`OPENCL.cl`** -- Top-level kernel entry point that assembles the above into a single compilable unit.

Kernel loading for exploration uses **opencltools** (`load_kernel_source`, `load_kernel_library`); runtime GPU assembly uses `kernel_loader.cpp` and `kernel_runners.cpp`, building on **opencltools** and **nmathopencl**.

### R-level infrastructure (`R/`)

40 files implementing the mixed-model formula-to-sampler pipeline (iid pfamily/prior/envelope/`rglmb`/`rlmb` R code has been removed; call `glmbayesCore::…` for that layer):

| File | Role |
|---|---|
| `model_setup.R` | lme4-style formula -> mixed-model design object (`model_setup()`) |
| `Prior_Setup_lmebayes.R` | Block~2 hyperprior calibration from reference `lmer` / `glmer` |
| `pfamily_list.R` / `pfamily_list_lmebayes_prior_setup.R` | S3 generic and `lmebayes_prior_setup` method -- Block~2 `pfamily` list from prior setup (embeds `glmbayesCore::dNormal()` / `dIndependent_Normal_Gamma()` objects) |
| `dGamma_list.R` / `dGamma_list_lmebayes_prior_setup.R` | Per-group Block~1 measurement-dispersion `glmbayesCore::dGamma()` pfamily list |
| `lme4_design_utilities.R` | Internal lme4 design chain (`get_lme4_components`, `extract_re_hyper_matrices`, …) |
| `rlmerb.R` / `rglmerb.R` | Matrix-level LMM / GLMM two-block samplers |
| `rLMM_reg.R` | Four Gaussian LMM replicate-chain routes (`rLMMNormal_reg*`, `rLMMindepNormalGamma_reg*`) plus dispatchers |
| `rGLMM_reg.R` / `rGLMM_sweep.R` | GLMM replicate-chain routes (known/estimated vcov) and the inner sweep-outer driver behind them |
| `two_block_rNormal_reg.R`, `two_block_batch_gibbs.R`, `two_block_pilot_cost.R`, `two_block_tau2_ref.R`, `two_block_measurement_prior.R`, `two_block_ergodicity.R`, `two_block_glmm_pilot_helpers.R`, `two_block_lmm_staged_sweep_outer.R` | Two-block Gibbs engine internals: Block~2 Normal regression, pilot-chain cost/TV calibration, dispersion reference tracking, ergodicity helpers |
| `two_block_sweep_history.R` / `plot_sweep_history_diag.R` | Sweep-history container, `print()` method, and cross-chain mean/SD diagnostic plot |
| `mixed_rmerb_helpers.R` | Internal helpers shared by `rlmerb()` / `rglmerb()` and **lmebayes** formula drivers |
| `build_mu_all.R` | Observation-level prior means for ICM / `simulate = FALSE` paths |
| `lmebayes_posterior_icm.R` | ICM posterior mean/mode (`lmerb_posterior_mean()`, `glmerb_posterior_mode()`) |
| `ing_prior_guard.R` | Truncation-window / ING dispersion prior guardrails |
| `rindepNormalGamma_reg_with_envelope.R` | `lmebayesCore`-specific envelope-based sampler (not a duplicate of any `glmbayesCore` export) |
| `lmebayes_reg_route_table.R` | Route dispatch table used by `.lmebayes_run_glmm_engine()` / `.lmebayes_run_lmm_engine()` |
| `glmmtmb_reference_helpers.R` | Optional `glmmTMB` reference-fit comparison helpers |

---

## Architecture: Two-Block Gibbs Sampling

`rlmerb()` / `rglmerb()` alternate between two conditional draws per outer sweep:

```
rlmerb(formula, pfamily_list, dispersion_ranef, data, n, ...)
       │
       ├─ 1. model_setup()              -- formula -> y / Z / groups / X_hyper
       ├─ 2. ICM posterior mean/mode      -- lmerb_posterior_mean() / build_mu_all()
       │      (always computed; used directly when simulate = FALSE)
       │
       └─ 3. simulate = TRUE: replicate-chain two-block Gibbs sweep
              │
              ├─ Block~1: group-level random effects b_j | gamma, tau^2
              │      (glmbayesCore::rglmb() / rlmb() per group, or the
              │       ING joint (b_j, sigma2_j) envelope path)
              │
              └─ Block~2: hyper means gamma_k | b, tau^2
                     via pfamily_list[[k]]$simfun -- Gaussian dNormal
                     components use glmbayesCore::rNormal_reg(); ING
                     dIndependent_Normal_Gamma components make a joint
                     (gamma_k, tau^2_k) draw via the envelope sampler
```

`rglmerb()` follows the same shape; non-Gaussian families route through `rGLMM_reg()` / `rGLMM_sweep()` (the sweep-outer replicate-chain driver, with an optional pilot stage governed by `gap_tol`) instead of the Gaussian `rLMM_reg` routes. Row-block engines (`block_rNormalReg()`, `block_rNormalGLM()`) apply the same two-block pattern independently within each row partition (used by **lmebayes**'s `lmbBlock()` / `glmbBlock()` via row-block priors, not by `rlmerb()`/`rglmerb()` directly).

All iid sampling inside a Block~1 or Block~2 draw (`rglmb()`, `rlmb()`, `rNormal_reg()`, the envelope machinery, and the `pfamily` prior constructors themselves) is **not** implemented in this package -- it is called via `glmbayesCore::…` at the R level, and resolved through the `glmbayesCore` namespace at the C++ level (see `src/package_ns.h`).

---

## Function overview

Symbols below are exported from **lmebayesCore** (`help(package = "lmebayesCore")`; 42 `export()` + 6 `S3method()` entries in `NAMESPACE`). End users typically load **lmebayes** instead, which re-exports the ones it needs directly.

### Model setup and priors

| Function | Role |
|----------|------|
| `model_setup()` | Parse an lme4-style formula into design matrices and variance components |
| `Prior_Setup_lmebayes()` | Calibrate Block~2 hyperpriors from a reference `lmer` / `glmer` fit |
| `pfamily_list()` | S3 generic; `pfamily_list.lmebayes_prior_setup()` builds Block~2 `pfamily` objects |
| `dGamma_list()` | S3 generic; `dGamma_list.lmebayes_prior_setup()` builds per-group Block~1 dispersion `pfamily` objects |
| `normalize_block()` | Row-block partition normalization (used by row-block engines and by **lmebayes**'s `lmbBlock()` / `glmbBlock()` / `Prior_SetupBlock()`) |
| `build_mu_all()` | Observation-level prior means when `simulate = FALSE` |

### Matrix-level two-block samplers

| Function | Role |
|----------|------|
| `rlmerb()` | Matrix-level Gaussian LMM two-block sampler (replicate chains) |
| `rglmerb()` | Matrix-level GLMM two-block sampler (`rGLMM_reg` routing for non-Gaussian; Gaussian delegates to `rLMM_reg` routes) |
| `lmerb_posterior_mean()` | Gaussian ICM fixef start when `simulate = FALSE` |
| `glmerb_posterior_mode()` | GLMM mode fixef start when `simulate = FALSE` |

### LMM / GLMM replicate-chain engines

| Function | Role |
|----------|------|
| `rLMMNormal_reg()`, `rLMMNormal_reg_known_vcov()`, `rLMMNormal_reg_estimated_vcov()` | Gaussian LMM routes with fixed `dispersion_ranef` |
| `rLMMindepNormalGamma_reg()`, `rLMMindepNormalGamma_reg_known_vcov()`, `rLMMindepNormalGamma_reg_estimated_vcov()` | Gaussian LMM routes with `dispersion_ranef = dGamma(...)` (ING) |
| `rGLMM_reg()`, `rGLMM_reg_known_vcov()`, `rGLMM_reg_estimated_vcov()` | Non-Gaussian GLMM replicate-chain dispatchers/routes |
| `rGLMM_sweep()` | Inner two-block sweep driver behind the `rGLMM_reg` and ING LMM routes |
| `rGLMM_Re_Draw()` | Re-draw helper for the GLMM sweep-outer engine |

### Row-block engines

| Function | Role |
|----------|------|
| `block_rNormalReg()`, `block_rNormalReg_update()` | Row-block Gaussian regression sampler and its incremental-update variant |
| `block_rNormalGLM()`, `block_rNormalGLM_update()` | Row-block GLM sampler and its incremental-update variant |

### Two-block internals and diagnostics

| Function | Role |
|----------|------|
| `two_block_rNormal_reg()` | Two-block Normal regression engine (Block~2 via `pfamily_list`) |
| `two_block_rate()`, `two_block_rate_from_pfamily_list()`, `two_block_tv_bound()`, `two_block_l_for_tv()` | Convergence-rate / total-variation calibration for inner Gibbs sweeps |
| `two_block_pilot_sampling_cost()`, `two_block_optimize_pilot_cost()`, `two_block_d0_pilot_start()`, `two_block_m_convergence_for_pilot_start()` | Pilot vs main chain cost optimization |
| `two_block_align_b_to_xhyper()` / `.two_block_align_b_to_xhyper` (+ `_cpp`) | Align Block~1 draws to `X_hyper` column order |
| `two_block_block2_one_chain()` / `.two_block_block2_one_chain` (+ `_cpp`) | Single-chain Block~2 update step |
| `plot_sweep_history_diag()` | Cross-chain mean/SD vs inner sweep for `two_block_sweep_history` |

### iid sampler retained in Core

| Function | Role |
|----------|------|
| `rindepNormalGamma_reg_with_envelope()` | Envelope-based joint coefficient/dispersion sampler specific to `lmebayesCore`'s calibration flow -- not a duplicate of any `glmbayesCore` export |

### S3 methods

`print.model_setup`, `print.lmebayes_prior_setup`, `print.two_block_sweep_history`, `print.two_block_rate`, `pfamily_list.lmebayes_prior_setup`, `dGamma_list.lmebayes_prior_setup`.

---

## Installation

**GitHub / R-Universe** (recommended for developers):

```r
install.packages("lmebayesCore",
                 repos = c("https://cloud.r-project.org",
                           "https://knygren.r-universe.dev"))
```

**From source** (required for OpenCL GPU support):

```r
# Ensure OpenCL development files are available on your system, then:
install.packages("lmebayesCore", type = "source",
                 repos = "https://knygren.r-universe.dev")
```

See [Chapter 16 -- Large models: GPU acceleration using OpenCL](https://knygren.r-universe.dev/articles/glmbayes/Chapter-16.html) for system-level setup instructions.

**Dependencies that must be installed first:**

```r
install.packages(c("Rcpp", "RcppArmadillo", "RcppParallel", "MASS", "Rdpack"))
install.packages("glmbayesCore",
                 repos = c("https://cloud.r-project.org",
                           "https://knygren.r-universe.dev"))
# opencltools and nmathopencl are available on R-Universe:
install.packages(c("opencltools", "nmathopencl"),
                 repos = "https://knygren.r-universe.dev")
```

---

## Extending lmebayesCore

### Mixed-model routes and pfamilies

Adding a new `pfamily` prior family (e.g. a new Block~1 or Block~2 conjugate/ING prior) is done in **glmbayesCore** -- see its `inst/ADDING_PFAMILY.md` -- since all `pfamily` constructors, `simfun` routing, and iid samplers now live there. `lmebayesCore` consumes those constructors through `pfamily_list()` / `dGamma_list()` and the two-block engines above.

### Block Gibbs ergodicity

`inst/BLOCK_GIBBS_ERGODICITY.md` documents the theoretical requirements for ergodicity when combining `glmbayesCore` iid simulation functions in the two-block Gibbs sampler implemented here.

### glmerb / two-block GLMM architecture

`inst/ARCHITECTURE_glmerb.md` maps the sweep-outer R driver (`rGLMM_sweep`), Block~1 / Block~2 call chains, formula-level wrappers in **lmebayes**, and the legacy C++ v5 path (removed) -- for maintainers working on `glmerb` parity or incremental C++ Block~2 ports.

### `R/` symbol inventory

Maintainers: exported API and internal helpers are listed in
[inst/R_FUNCTION_INVENTORY.md](inst/R_FUNCTION_INVENTORY.md)
(`R_EXPORTED_AND_DOCUMENTED.md`, `R_INTERNAL_HELPERS.md`). **Note:** these
inventory files predate the `glmbayesCore` deduplication (Stages 0-1c) and
still describe the pre-fork symbol set; treat them as historical context
pending a refresh, not as an accurate current-state reference.

---

## Key References

- Nygren, K.N. and Nygren, L.M. (2006). Likelihood Subgradient Densities. *Journal of the American Statistical Association*, 101(475), 1144-1156. -- The accept-reject envelope method underlying **glmbayesCore**'s non-Gaussian iid samplers, called from Block~1 / Block~2 here.
- Lindley, D.V. and Smith, A.F.M. (1972). Bayes estimates for the linear model. *Journal of the Royal Statistical Society B*, 34, 1-41. -- Conjugate Normal-Gamma foundations.
- Gelman, A. et al. (2013). *Bayesian Data Analysis*, 3rd ed. -- Reference for prior specifications and dispersion modeling.

A complete bibliography is in `inst/REFERENCES.bib`.

---

## Future plans

- **Stage 3 -- C++ deduplication:** `lmebayesCore`'s `src/` still compiles its own copy of the full iid envelope/sampler engine (the tables under [What Is Inside lmebayesCore](#what-is-inside-lmebayesCore) above), even though the R layer no longer exposes any of it directly -- the mixed-model `.cpp` files (`twoBlockGibbs.cpp` and friends) call into that object code within the same DLL. A future pass should audit which of the ~34 shared `.cpp` files are safe to exclude from this package's build (via `src/Makevars` or deletion) once verified unused outside the mixed-model paths that still need them (e.g. envelope/fam/sim code used by `block_rIndepNormalGammaReg.cpp`, `rNormalRegBlocks.cpp`, `two_block_block1.cpp`, and Block~1 in `twoBlockGibbs.cpp`).
- **Sweep-outer drivers and `sweep_history` on all two-block paths:** Mixed-model
  sampling should eventually use a **sweep-outer** loop (all chains complete inner
  sweep `m`, then `m+1`, …) on every route, for consistency with `rGLMM_sweep()` /
  `rGLMM_reg_*`. Each stored draw should attach **`sweep_history`** (class
  `two_block_sweep_history`) so `print()` and `plot_sweep_history_diag()` can
  diagnose inner-Gibbs convergence. Today, sweep-outer R drivers and ING pilot/main
  paths already capture history; **gaps remain** (e.g. Gaussian `lmerb()` with fixed
  σ² and fixed Block~2 τ² still uses the v2 C++ chain-outer driver, which does not
  export per-sweep cross-chain stats). See `inst/ARCHITECTURE_glmerb.md` and
  `inst/PLAN_block2_cpp_migration.md`.
- **C++ inner-chain loops and within-block parallel sampling:** The per-sweep
  **inner chain** loops (Block~1 + Block~2 updates across replicate chains) should
  migrate from R orchestration (`rGLMM_sweep()`, batch helpers in
  `two_block_batch_gibbs.R`) into **`src/*.cpp`** drivers (alongside / replacing
  v2 and v5), matching the v5 sweep-outer layout. Parallelism should be
  **within-block, across chains** at fixed inner sweep `m` (not parallel inner
  sweeps): **Block~1** random-effect updates over replicate chains first
  (**higher priority**); **Block~2** fixed-effect / hyperparameter updates over
  chains in parallel where safe (**ideal follow-on**). Use native threading (e.g.
  `RcppParallel`) so large `n` does not pay full R-loop overhead.

---

## License

GPL-2. See the `LICENSE` file and `inst/COPYRIGHTS` for attribution of incorporated R Mathlib sources.
