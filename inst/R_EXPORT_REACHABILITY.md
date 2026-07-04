# Export reachability — dead / inactive analysis

Which **glmbayesCore** `@export` symbols are on a live call path from
end-user packages, vs candidates for `@noRd` / unexport review.

**Companions:** [R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md)
(overlap matrix), [R_CORE_ONLY_EXPORTS.md](R_CORE_ONLY_EXPORTS.md) (by function type).

**Last reviewed:** 2026-07-04.

---

## Reachability rules

A Core export is **live** if it is reached (directly or transitively in Core `R/`)
from either:

1. **glmbayes planned retain re-exports** (13 symbols — stay on
   `library(glmbayes)` after migration), or
2. **lmebayes** exported drivers and Core symbols re-exported / imported there.

**Phase-out** symbols (simfunctions, envelopes, `*_ct`, …) are **live** when
**`rglmb()` / `rlmb()`** reach them via `pfamily$simfun`, even though
**glmbayes** will drop them from `NAMESPACE` and callers will use
`glmbayesCore::`.

**Not live:** no caller on the graph from (1) or (2). User composition
(e.g. `multi_prior_setup()` then `multi_rlmb()`) counts once **`multi_rlmb()`**
is a planned retain re-export.

**Out of scope:** `data-raw/`, `inst/examples/`, demos, and direct
`glmbayesCore::` from maintainer scripts (may still justify keeping exports).

---

## Entry points

### glmbayes — planned retain re-exports (13)

| Function | Typical Core reach |
|----------|-------------------|
| `Prior_Setup()` | `compute_gaussian_prior()` |
| `Prior_Check()` | `glm()` only (no sim stack) |
| `pfamily()`, `dNormal()`, `dGamma()`, `dBeta()`, `dNormal_Gamma()`, `dIndependent_Normal_Gamma()` | Embedded `simfun` pointers |
| `multi_prior_setup()` | `Prior_Setup()` per `cbind` column |
| `multi_rlmb()` | `rlmb()` per column (pairs with `multi_prior_setup()`) |
| `rglmb()`, `rlmb()` | `pfamily$simfun` → simfunction / envelope pipeline |
| `diagnose_glmbayes()` | `glmbayesCore_has_opencl()` |

**Note:** `multi_rlmb()` is not yet in **glmbayes** `NAMESPACE`; policy is to
add it as a retain re-export. Same help topic exports
`multi_rNormalGamma_reg()` and `multi_rindepNormalGamma_reg()` — live via that
API (each calls phase-out `rNormalGamma_reg()` / `rindepNormalGamma_reg()` per
column).

### lmebayes — additional entry points

| Function | Role |
|----------|------|
| `lmerb()`, `glmerb()` | Formula mixed-model drivers |
| `rlmerb()`, `rglmerb()` | Re-exported matrix samplers |
| `model_setup()`, `Prior_Setup_lmebayes()`, `pfamily_list()` | Design + Block~2 priors |
| `Prior_SetupBlock()`, `lmbBlock()`, `glmbBlock()` | Row blocks (`normalize_block`; block fits via **glmbayes** `lmb`/`glmb` → planned Core `rlmb`/`rglmb`) |
| `plot_sweep_history_diag()` | Sweep-history plots |
| Re-exported retain symbols | Same as glmbayes retain subset |

**Direct Core calls (not re-exported from lmebayes):** `build_mu_all()`,
`lmerb_posterior_mean()`, `glmerb_posterior_mode()`, `normalize_block()`.

---

## Live paths (summary)

### From planned glmbayes retain

```
Prior_Setup() → compute_gaussian_prior()
multi_prior_setup() → Prior_Setup() (× columns)

multi_rlmb() → rlmb() (× columns) → pfamily$simfun → …
multi_rNormalGamma_reg() → rNormalGamma_reg() (× columns)
multi_rindepNormalGamma_reg() → rindepNormalGamma_reg() (× columns)

rglmb() / rlmb() → pfamily$simfun →
  rNormal_reg(), rNormalGamma_reg(), rindepNormalGamma_reg(),
  rGamma_reg(), rGamma_Conjugate_reg(), rBeta_reg(),
  rNormalGLM_std(), rIndepNormalGammaReg_std(),
  glmbfamfunc(), glmb.wfit(), rNormal_reg.wfit(),
  glmb_Standardize_Model(), EnvelopeSize … EnvelopeOrchestrator,
  pnorm_ct(), rnorm_ct(), pinvgamma_ct(), qinvgamma_ct(),
  rinvgamma_ct(), rgamma_ct(), …

diagnose_glmbayes() → glmbayesCore_has_opencl()
```

### From lmebayes (mixed models)

```
lmerb() / glmerb()
  → rlmerb() / rglmerb()
    → rLMMNormal_reg() | rLMMNormal_reg_estimated_vcov()
      | rLMMindepNormalGamma_reg() | rGLMM()
        → two_block_rNormal_reg(), rGLMM_sweep(), rGLMM_Re_Draw()
        → block_rNormalReg() / block_rNormalGLM() (Block~1)
        → rglmb() (Block~2 ING / hyper draws)
        → two_block_rate_from_pfamily_list(), two_block_rate(), two_block_l_for_tv(),
           two_block_tv_bound(), two_block_mode_weights(),
           two_block_optimize_pilot_cost(), two_block_pilot_* , …
        → two_block_block2_one_chain_cpp() [default];
           two_block_block2_one_chain(), two_block_align_b_to_xhyper() [R fallback]
```

---

## Flagged: dead or inactive exports

### High confidence — not on retain or lmebayes function graph

| Export | File (approx.) | Why inactive |
|--------|----------------|--------------|
| `multi_rNormal_reg()` | `multi_rNormal_reg.R` | Parallel matrix API; not called by `multi_rlmb()` or retain/lmebayes drivers. |
| `block_rNormalReg_update()` | `rNormalReg_reg_block_update.R` | Not on `rGLMM_sweep` / `lmerb` chain; examples / roxygen only. |
| `block_rNormalGLM_update()` | `rNormalGLM_reg_block_update.R` | Same. |
| `simfunction()` | `simfunction.R` | Introspection generic; runtime uses `pfamily$simfun` directly. |

**Note:** `block_rNormalReg()` and `block_rNormalGLM()` are **live** (called from
`.two_block_block1_draw_block()` on non-Gaussian **glmerb** paths). Only the
`*_update` wrappers are flagged.

### Medium confidence — default-off or niche

| Export | Status |
|--------|--------|
| `two_block_block2_one_chain()`, `two_block_align_b_to_xhyper()` (R) | Live only when `rGLMM_sweep(..., use_cpp_block2 = FALSE)`. Default C++ path. |
| `.two_block_align_b_to_xhyper`, `.two_block_block2_one_chain` | Namespace aliases; same. |

---

## Explicitly not dead (common false positives)

| Export | Why live |
|--------|----------|
| Phase-out sim / envelope / `*_ct` (30 symbols) | Reachable from planned **`rglmb()` / `rlmb()`** |
| `compute_gaussian_prior()` | **`Prior_Setup()`** (retain) |
| `glmbayesCore_has_opencl()` | **`diagnose_glmbayes()`** (retain) |
| `multi_rlmb()`, `multi_rNormalGamma_reg()`, `multi_rindepNormalGamma_reg()` | Planned **glmbayes** retain / co-export multi-response API |
| Core `rlmb()`, `rglmb()`, `multi_prior_setup()` | Planned retain |
| `rGLMM()`, `rGLMM_sweep()`, `rGLMM_Re_Draw()`, `rLMMNormal_reg*`, `two_block_rNormal_reg()`, TV/pilot helpers | **lmebayes** `lmerb()` / `glmerb()` chain |
| `rLMMNormal_reg_known_vcov()` | Routed from `rLMMNormal_reg()` |

---

## Action shortlist (maintainers)

Candidates for `@noRd` / unexport after `data-raw` audit:

1. `multi_rNormal_reg` (not part of `multi_rlmb()` chain)
2. `block_rNormalReg_update` + `block_rNormalGLM_update`
3. `simfunction()` generic

Do **not** unexport on reachability grounds alone: phase-out sim/envelope
exports, mixed-model stack, lmebayes direct imports (`build_mu_all`, …), or
planned retain **`multi_rlmb()`** and its co-export siblings.

---

## Review checklist

| Priority | Item |
|----------|------|
| 1 | Re-run this graph after **glmbayes** `Imports: glmbayesCore` lands and **`multi_rlmb()`** is re-exported. |
| 2 | Before dropping legacy exports, grep `data-raw/` and vignettes. |
| 3 | Keep R Block~2 fallbacks until C++ path is universal and debug scripts migrated. |
