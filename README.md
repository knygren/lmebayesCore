# glmbayesCore

![GitHub release (latest by date)](https://img.shields.io/github/v/release/knygren/glmbayesCore?label=version)
![License: GPL-2](https://img.shields.io/badge/license-GPL--2-blue.svg)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/knygren/glmbayesCore/R-CMD-check.yaml?label=R%20CMD%20Check)

**glmbayesCore** is the compiled sampling engine that powers the glmbayes ecosystem. It holds the C++/OpenCL samplers, the family-function infrastructure, and the R-level prior and simulation interfaces that downstream packages depend on. End users should install [glmbayes](https://github.com/knygren/glmbayes) rather than this package directly.

The relationship to the broader ecosystem parallels how `StanHeaders` / `rstan` serve as the compiled backbone for `rstanarm`: glmbayesCore is the infrastructure layer; `glmbayes` (and the in-development `lmebayes`) are the user-facing packages built on top of it.

---

## Package Ecosystem

```
                ┌─────────────────────────────────────────┐
                │           End-user packages             │
                │   glmbayes  ·  lmebayes  ·  (others)   │
                └──────────────────┬──────────────────────┘
                                   │ Imports / LinkingTo
                ┌──────────────────▼──────────────────────┐
                │              glmbayesCore               │
                │  C++ samplers · OpenCL kernels          │
                │  pfamily · simfunctions · rglmb/rlmb   │
                └──────────────────┬──────────────────────┘
                                   │ Imports
                ┌──────────────────▼──────────────────────┐
                │   opencltools  ·  nmathopencl            │
                │   Rcpp · RcppArmadillo · RcppParallel   │
                └─────────────────────────────────────────┘
```

**glmbayes** adds the formula interface (`glmb()`, `lmb()`), MCMC diagnostics, and the full suite of S3 methods that mirror base-R's `lm()` / `glm()`.

**lmebayes** (in development) extends the engine to linear mixed-effects models, interfacing with glmbayesCore at both the C++ and R levels.

---

## What Is Inside glmbayesCore

### C++ sampling engine (`src/`)

The core is organized under the `glmbayes::` namespace, partitioned into sub-namespaces:

| Sub-namespace | Key files | Role |
|---|---|---|
| `glmbayes::fam` | `famfuncs.h`, `famfuncs_*.cpp` | Negative log-posterior (`f2`) and gradient (`f3`) for gaussian, poisson, binomial, Gamma |
| `glmbayes::env` | `EnvelopeBuild*.cpp`, `EnvelopeEval.cpp`, `EnvelopeSort.cpp`, `EnvelopeSize.cpp`, `Set_Grid.cpp`, `Set_LogP.cpp` | Piecewise-exponential envelope construction (Nygren & Nygren, 2006) |
| `glmbayes::sim` | `rNormalGLM.cpp`, `rIndepNormalGammaReg.cpp`, `rNormalGammaReg.cpp`, `rNormalReg.cpp`, `rGammaGamma.cpp`, `rGammaGaussian.cpp` | Posterior samplers |
| `glmbayes::rng` | `rng_utils.cpp` | Thread-safe RNG wrappers for parallel sampling |
| `glmbayes::progress` | `progress_utils.cpp` | Optional progress bar support |

Export wrappers in `export_wrappers.cpp` and `kernel_wrappers.cpp` expose selected entry points to R via Rcpp.

### OpenCL kernels (`inst/cl/`)

For systems with an OpenCL-capable device, envelope construction can be offloaded to the GPU. The `inst/cl/` tree contains:

- **`src/f2_f3_*.cl`** — GPU ports of the `f2` (negative log-posterior) and `f3` (gradient) functions for each family/link combination (gaussian, poisson, binomial logit/probit/cloglog, Gamma log).
- **`nmath/`** — OpenCL port of the R Mathlib probability functions (`dnorm`, `dgamma`, `dbinom`, `dpois`, `lgamma`, etc.) needed by the kernels.
- **`libR_shims/`, `R_ext_*/`, `R_shims/`** — Shim headers that make the nmath kernels compile cleanly under OpenCL's C99-based dialect.
- **`OPENCL.cl`** — Top-level kernel entry point that assembles the above into a single compilable unit.

Kernel loading for exploration uses **opencltools** (`load_kernel_source`, `load_kernel_library` with `package = "glmbayesCore"`); runtime GPU assembly uses `kernel_loader.cpp` and `kernel_runners.cpp`, building on **opencltools** and **nmathopencl**.

### R-level infrastructure (`R/`)

| File | Role |
|---|---|
| `pfamily.R` | Prior-family constructors (`dNormal`, `dNormal_Gamma`, `dIndependent_Normal_Gamma`, `dGamma`, `dBeta`) and the `pfamily()` generic |
| `prior.R` | `Prior_Setup()`, `Prior_Check()`, and helper utilities for computing default hyperparameters |
| `simfunction.R` | Low-level simulation functions (`rNormal_reg`, `rNormalGamma_reg`, `rindepNormalGamma_reg`, `rGamma_reg`, `rGamma_Conjugate_reg`, `rBeta_reg`) and the `simfunction()` introspection generic |
| `simulationpipeline.R` | `glmbfamfunc()` (R closure bundle for f1–f4 and f7), pipeline documentation, and fit helpers |
| `rglmb.R` / `rlmb.R` | Matrix-input samplers — the primary R-level interface consumed by downstream packages |
| `pfamily_list.R` / `pfamily_list_lmebayes_prior_setup.R` | S3 generic and `lmebayes_prior_setup` method — Block~2 `pfamily` list from prior setup |
| `model_setup.R` | lme4-style formula → mixed-model design object (`model_setup()`) |
| `Prior_Setup_lmebayes.R` | Block~2 hyperprior calibration from reference `lmer` / `glmer` |
| `lme4_design_utilities.R` | Internal lme4 design chain (`get_lme4_components`, `extract_re_hyper_matrices`, …) |
| `rlmerb.R` / `rglmerb.R` | Matrix-level LMM / GLMM two-block samplers |
| `mixed_rmerb_helpers.R` | Internal helpers for `rlmerb()` / `rglmerb()` and **lmebayes** formula drivers |
| `plot_sweep_history_diag.R` | Cross-chain mean/SD plots for `two_block_sweep_history` (pilot/main stages) |
| `two_block_sweep_history.R` | Sweep-history container and `print.two_block_sweep_history()` |
| `envelopeorchestrator.R` | R orchestration of multi-step envelope building and optional GPU dispatch |
| `compute_gaussian_prior.R` | Gaussian-specific prior calibration utilities |

---

## Architecture: How pfamilies Route to Simulation Functions

A `pfamily` object is a self-contained prior specification. Every constructor bundles the hyperparameters into a `prior_list` **and** embeds a `simfun` function pointer that knows how to sample the corresponding posterior. When `rglmb()` (or any downstream modelling function) draws samples, it simply calls `pfamily$simfun(y, x, prior_list, family, ...)` — there is no internal `switch` on prior type.

```
rglmb(y, x, pfamily = dNormal(...), family = poisson())
          │
          └─► pfamily$simfun  ──►  rNormal_reg()
                                       │
                              family == gaussian?
                              ├── Yes ──► conjugate multivariate normal draw
                              └── No  ──► envelope sampling (Nygren & Nygren, 2006)
                                              │
                                              └──► rNormalGLM (C++)
```

The full routing table for the implemented pfamilies:

| pfamily constructor | Embedded `simfun` | Posterior path |
|---|---|---|
| `dNormal()` | `rNormal_reg()` | Conjugate MVN draw (Gaussian); subgradient envelope sampling (all other families) |
| `dNormal_Gamma()` | `rNormalGamma_reg()` | Conjugate Normal-Gamma draw (Gaussian only) |
| `dIndependent_Normal_Gamma()` | `rindepNormalGamma_reg()` | Joint coefficient + dispersion envelope (Gaussian; non-conjugate) |
| `dGamma(Inv_Dispersion = TRUE)` | `rGamma_reg()` | Gamma prior on inverse dispersion; accept-reject or conjugate draw |
| `dGamma(Inv_Dispersion = FALSE)` | `rGamma_Conjugate_reg()` | Conjugate Gamma–Poisson or Gamma–Gamma (intercept-only, identity link) |
| `dBeta()` | `rBeta_reg()` | Conjugate Beta–Binomial (intercept-only, identity link) |

`Prior_Setup()` fits an auxiliary GLM and returns calibrated hyperparameters (`mu`, `Sigma`, `shape`, `rate`, etc.) on the same scale as the design matrix. Its output slots into any pfamily constructor directly.

---

## Architecture: How Simulation Functions Route to C++ Samplers

The R-level simulation functions are thin orchestration wrappers. Their main jobs are: validate and pre-process inputs, select the correct C++ entry point, and post-process the returned list into a classed object. Several simulation functions dispatch to more than one C++ sampler depending on the family or the model state.

### `rNormal_reg()`

```
rNormal_reg(y, x, prior_list, family, ...)
       │
  family$family == "gaussian"?
  ├── Yes ──► direct MVN draw via backsolve / Cholesky   (pure R / RcppArmadillo)
  └── No  ──► EnvelopeOrchestrator (R)
                   ├── EnvelopeBuild (C++)          [construct piecewise envelope]
                   │       uses famfuncs_*.cpp for f2/f3 per family
                   └── rNormalGLM (C++)             [accept-reject sampling]
                               ├── RcppParallel worker threads (TBB)
                               └── optional OpenCL path for envelope build
```

### `rindepNormalGamma_reg()`

This simulation function handles a *joint* prior over regression coefficients **and** dispersion. The dispersion enters the envelope through a grid over plausible dispersion values, so the sampler loops over that grid inside a single C++ call:

```
rindepNormalGamma_reg(y, x, prior_list, ...)
       │
       └──► rIndepNormalGammaReg (C++)
                   ├── for each dispersion grid point:
                   │       EnvelopeBuild_Ind_Normal_Gamma (C++)
                   │           [conditional coefficient envelope at fixed dispersion]
                   └── joint accept-reject over (beta, dispersion) pairs
```

### `rGamma_reg()`

```
rGamma_reg(y, x, prior_list, family, ...)
       │
  family$family == "gaussian"?
  ├── Yes ──► rGammaGaussian (C++)   [Gamma draw on precision of Gaussian model]
  └── No  ──► rGammaGamma (C++)     [Gamma draw on dispersion of Gamma(log) model]
```

---

## Architecture: How `rglmb()` Orchestrates a Draw

`rglmb()` is the canonical matrix-input orchestrator for GLM posterior draws. It does not sample directly; instead it validates the `family × pfamily` combination and delegates all sampling work to the simulation function that the `pfamily` object carries. In the **glmbayes** package, `glmb()` and `lmb()` wrap `rglmb()` with formula parsing and model-frame construction; downstream packages can build similar orchestrators on the same pattern.

The sequence for a single `rglmb()` call is:

```
rglmb(y, x, family = poisson(), pfamily = dNormal(mu, Sigma), n = 1000)
  │
  ├─ 1. Resolve family
  │       normalise string / function → family object
  │       special case: Poisson + dGamma(Inv_Dispersion=TRUE) → coerce to conjugate rate prior
  │
  ├─ 2. Unpack pfamily
  │       okfamilies  ← which family$family values are allowed
  │       plinks      ← function: family → allowed link strings
  │       prior_list  ← hyperparameters (mu, Sigma, shape/rate, …)
  │       simfun      ← pointer to the backend sampler (e.g. rNormal_reg)
  │
  ├─ 3. Validate the combination
  │       family$family ∈ okfamilies?   → error if not
  │       family$link   ∈ plinks(family)?  → error if not
  │       (extra geometry checks for scalar-only conjugate priors)
  │
  ├─ 4. Call the sampler
  │       outlist ← simfun(n, y, x, prior_list, family, offset, weights,
  │                         Gridtype, n_envopt, use_parallel, use_opencl, …)
  │
  └─ 5. Post-process and return
          overwrite call slot with match.call()
          re-attach pfamily and simfun_args (for traceability)
          set coefficient and coef.mode names from colnames(x)
          return object of class c("rglmb", "glmb", "glm", "lm")
```

The key design point is step 4: **`rglmb()` contains no `switch` on prior type**. The right sampler was bound to `simfun` at the moment the user called `dNormal()` / `dGamma()` / etc., so `rglmb()` simply calls whatever function is sitting there. Adding a new prior family therefore requires no changes to `rglmb()` itself — only a new pfamily constructor and a new simulation function.

The interaction with the two architecture sections above is:

- **pfamily routing table** — determines which `simfun` is embedded (step 2).
- **Simulation function → C++ routing** — what happens inside `simfun(...)` (step 4).

`rlmb()` follows the same pattern but restricts `okfamilies` to `"gaussian"` and skips the `glmbfamfunc` step, since the Gaussian posterior is always conjugate or near-conjugate.

**Extensibility note.** The orchestrator pattern is intentionally generic: validate a model specification, unpack a routing object, call the embedded function, post-process. A mixed-effects package such as **lmebayes** implements formula drivers (`lmerb()`, `glmerb()`) on top of `rlmerb()` / `rglmerb()` and re-exports the mixed-model setup symbols below. The only requirement at each Gibbs block is a compatible `prior_list` for the relevant `simfun`.

---

## Function overview

Symbols below are exported from **glmbayesCore** (`help(package = "glmbayesCore")`). End users typically load **glmbayes** or **lmebayes** instead; those packages re-export subsets of this API.

**Maintainers:** full export and helper inventories live in
[inst/R_FUNCTION_INVENTORY.md](inst/R_FUNCTION_INVENTORY.md)
([exports / overlap matrix](inst/R_EXPORTED_AND_DOCUMENTED.md),
[Core-only by function type](inst/R_CORE_ONLY_EXPORTS.md),
[export reachability](inst/R_EXPORT_REACHABILITY.md),
[internal helpers](inst/R_INTERNAL_HELPERS.md)).

### Shared with **glmbayes** (iid GLM / LM)

**glmbayes** currently re-exports 42 symbols from **glmbayesCore**. Maintainer
policy: keep the user-facing prior/sampler API on `library(glmbayes)`; phase
low-level simulation and envelope exports to **glmbayesCore**-only (see
[inst/R_EXPORTED_AND_DOCUMENTED.md](inst/R_EXPORTED_AND_DOCUMENTED.md)).

#### Retain as **glmbayes** re-exports

| Function | Role |
|----------|------|
| `Prior_Setup()`, `Prior_Check()` | Default prior calibration and prior predictive checks |
| `pfamily()`, `dNormal()`, `dNormal_Gamma()`, `dIndependent_Normal_Gamma()`, `dGamma()`, `dBeta()` | Prior-family constructors |
| `multi_prior_setup()` | Multi-response Gaussian prior setup (`cbind` LHS; calls `Prior_Setup()` per column) |
| `multi_rlmb()` | Multi-response LM sampler (`cbind` LHS; `rlmb()` per column). Planned re-export — not yet in **glmbayes** `NAMESPACE`. |
| `rglmb()`, `rlmb()` | Matrix-level Bayesian GLM / LM samplers (`glmb()` / `lmb()` backends) |
| `diagnose_glmbayes()` | OpenCL / GPU diagnostic report |

#### Phase out of **glmbayes** (stay in **glmbayesCore**)

| Function | Role |
|----------|------|
| `compute_gaussian_prior()` | Internal Gaussian calibration used only inside `Prior_Setup()` |
| `simfunction()`, `glmbfamfunc()` | Simulation registry and GLM family pipeline helpers |
| `rNormal_reg()`, `rNormalGamma_reg()`, `rindepNormalGamma_reg()`, `rGamma_reg()`, `rBeta_reg()`, … | Low-level `simfunction` samplers |
| `rNormalGLM_std()`, `rIndepNormalGammaReg_std()`, `glmb.wfit()`, `glmb_Standardize_Model()` | Standardized envelope path and fitter hooks |
| `EnvelopeBuild()`, `EnvelopeOrchestrator()`, `EnvelopeSize()`, … | Accept–reject envelope machinery |
| `pnorm_ct()`, `rnorm_ct()`, `pinvgamma_ct()`, `rgamma_ct()`, … | Truncated-distribution C++ callbacks |

See the **glmbayes** README and vignettes for the formula interface (`glmb()`, `lmb()`) and S3 methods built on the retained exports.

### Mixed-model setup and sampling (also re-exported by **lmebayes**)

| Function | Role |
|----------|------|
| `model_setup()` | Parse an lme4-style formula into design matrices and variance components |
| `Prior_Setup_lmebayes()` | Calibrate Block~2 hyperpriors from a reference `lmer` / `glmer` fit |
| `pfamily_list()` | S3 generic; `pfamily_list.lmebayes_prior_setup()` builds Block~2 `pfamily` objects |
| `rlmerb()` | Matrix-level Gaussian LMM two-block sampler (replicate chains) |
| `rglmerb()` | Matrix-level GLMM two-block sampler (`rLMMNormal_reg` / `rGLMM` routing) |
| `plot_sweep_history_diag()` | Cross-chain mean/SD vs inner sweep for `two_block_sweep_history` |

Typical **lmebayes** workflow: `model_setup()` → `Prior_Setup_lmebayes()` → `pfamily_list(ps)` → `lmerb()` / `glmerb()`.

S3 helpers: `print.model_setup`, `print.lmebayes_prior_setup`, `print.two_block_sweep_history`.

### **lmebayes** direct Core calls (`importFrom` or `glmbayesCore::` — must stay exported)

| Function | **lmebayes** callers | Role |
|----------|----------------------|------|
| `build_mu_all()` | `lmerb()`, `glmerb()` | Observation-level prior means when `simulate = FALSE` |
| `lmerb_posterior_mean()` | `lmerb()` | Gaussian ICM fixef start when `simulate = FALSE` |
| `glmerb_posterior_mode()` | `glmerb()` | GLMM mode fixef start when `simulate = FALSE` |
| `normalize_block()` | `lmbBlock()`, `glmbBlock()`, `Prior_SetupBlock()` | Row-block partition normalization |

### Two-block engines — indirect from **lmebayes** (export optional for **lmebayes**)

| Function | Role |
|----------|------|
| `two_block_rNormal_reg_v2()` | Two-block Normal regression engine (Block~2 via `pfamily_list`) |
| `two_block_rate_v2()`, `two_block_mode_weights()`, `two_block_tv_bound()` | TV / rate calibration for inner Gibbs sweeps |
| `two_block_optimize_pilot_cost()` | Pilot vs main chain cost optimization |
| `rGLMM()`, `rGLMM_sweep()` | GLMM sweep-outer driver; **lmebayes** reaches via `glmerb()` → `rglmerb()` only |
| `rLMMNormal_reg()`, `rLMMindepNormalGamma_reg()`, `rLMMNormal_reg_estimated_vcov()` | Gaussian LMM routers; **lmebayes** reaches via `lmerb()` / `glmerb()` → `rlmerb()` / `rglmerb()` only |
| `block_rNormalReg()`, `block_rNormalGLM()` | Row-block samplers for BY-style splits |

These are listed under **glmbayesCore-only exports** in `inst/R_EXPORTED_AND_DOCUMENTED.md` (indirect from **lmebayes** subsection). Export is optional for **lmebayes** — Core routes inside `rlmerb()` / `rglmerb()`.

### Internal helpers

Undocumented `@noRd` symbols (mixed-model glue, lme4 design chain, two-block staging, envelope internals) are listed in [inst/R_INTERNAL_HELPERS.md](inst/R_INTERNAL_HELPERS.md). **lmebayes** resolves a subset via `glmbayesCore:::` / `importFrom`.

---

## Developer Interface Levels

Downstream packages and developers can interface with glmbayesCore at several levels:

### Level 1 — C++ (via `LinkingTo`)

Add `glmbayesCore` to `LinkingTo:` in your `DESCRIPTION` and include the exported headers:

```cpp
#include "glmbayesCore/famfuncs.h"      // f2 / f3 for all families
#include "glmbayesCore/Envelopefuncs.h" // envelope build / eval routines
#include "glmbayesCore/simfuncs.h"      // glmb_Standardize_Model and friends
#include "glmbayesCore/R_interface.h"   // R ↔ C++ data conversion helpers
```

This is the lowest-level entry point and gives full access to the compiled samplers and envelope machinery without going through R.

### Level 2 — R simulation functions

Call the simulation functions directly from R, bypassing the formula interface entirely. This is well-suited to Gibbs samplers and other workflows where you hold some parameters fixed and update others:

```r
library(glmbayesCore)

# Full joint draw: coefficients + dispersion under independent Normal-Gamma prior
fit <- rindepNormalGamma_reg(
  y          = y,
  x          = X,
  n          = 2000,
  prior_list = dIndependent_Normal_Gamma(mu, Sigma, shape, rate)$prior_list,
  family     = gaussian()
)
```

### Level 3 — `rglmb()` / `rlmb()` with pfamily objects

The matrix-input samplers `rglmb()` and `rlmb()` are the canonical R interface. They accept any pfamily object and handle all routing automatically:

```r
ps  <- Prior_Setup(y, X, family = poisson())
fit <- rglmb(y = y, x = X, n = 1000,
             pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
             family  = poisson())
```

Downstream packages such as `glmbayes` wrap this level with formula parsing, model-frame construction, and the full S3 method infrastructure.

---

## Installation

**GitHub / R-Universe** (recommended for developers):

```r
install.packages("glmbayesCore",
                 repos = c("https://cloud.r-project.org",
                           "https://knygren.r-universe.dev"))
```

**From source** (required for OpenCL GPU support):

```r
# Ensure OpenCL development files are available on your system, then:
install.packages("glmbayesCore", type = "source",
                 repos = "https://knygren.r-universe.dev")
```

See [Chapter 16 — Large models: GPU acceleration using OpenCL](https://knygren.r-universe.dev/articles/glmbayes/Chapter-16.html) for system-level setup instructions.

**Dependencies that must be installed first:**

```r
install.packages(c("Rcpp", "RcppArmadillo", "RcppParallel", "MASS", "Rdpack"))
# opencltools and nmathopencl are available on R-Universe:
install.packages(c("opencltools", "nmathopencl"),
                 repos = "https://knygren.r-universe.dev")
```

---

## Extending glmbayesCore

### Adding a new pfamily

The file `inst/ADDING_PFAMILY.md` contains a step-by-step guide. In summary:

1. Write a constructor in `pfamily.R` that builds `prior_list` and sets `simfun`.
2. Implement or reuse a simulation function in `simfunction.R`.
3. If a new C++ sampler is needed, add it under `src/`, register it via `RcppExports`, and expose it through `export_wrappers.cpp`.
4. For GPU support, add the corresponding `f2`/`f3` OpenCL kernel under `inst/cl/src/` and register it in the assembly path in `kernel_loader.cpp`.

### Block Gibbs ergodicity

`inst/BLOCK_GIBBS_ERGODICITY.md` documents the theoretical requirements for ergodicity when combining glmbayesCore simulation functions in a block Gibbs sampler — relevant for packages like `lmebayes` that cycle between coefficient and variance-component updates.

### glmerb / two-block GLMM architecture

`inst/ARCHITECTURE_glmerb.md` maps the sweep-outer R driver (`rGLMM_sweep`), Block 1 / Block 2 call chains, formula-level wrappers in **lmebayes**, and the legacy C++ v5 path — for maintainers working on `glmerb` parity or incremental C++ Block 2 ports.

### `R/` symbol inventory

Maintainers: exported API and internal helpers are listed in
[inst/R_FUNCTION_INVENTORY.md](inst/R_FUNCTION_INVENTORY.md)
(`R_EXPORTED_AND_DOCUMENTED.md`, `R_INTERNAL_HELPERS.md`).

---

## Key References

- Nygren, K.N. and Nygren, L.M. (2006). Likelihood Subgradient Densities. *Journal of the American Statistical Association*, 101(475), 1144–1156. — The accept-reject envelope method at the heart of the non-Gaussian samplers.
- Lindley, D.V. and Smith, A.F.M. (1972). Bayes estimates for the linear model. *Journal of the Royal Statistical Society B*, 34, 1–41. — Conjugate Normal-Gamma foundations.
- Gelman, A. et al. (2013). *Bayesian Data Analysis*, 3rd ed. — Reference for prior specifications and dispersion modeling.

A complete bibliography is in `inst/REFERENCES.bib`.

---

## License

GPL-2. See the `LICENSE` file and `inst/COPYRIGHTS` for attribution of incorporated R Mathlib sources.
