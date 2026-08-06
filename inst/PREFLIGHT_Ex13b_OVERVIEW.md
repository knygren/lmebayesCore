# Ex_13b preflight: four functions before the sampler

`demo("Ex_13b_rLMM_estimated_dispersion_known_vcov_BigWordClub_PartVI", package = "lmebayesCore")` builds design and priors with four exported entry points, then calls `rLMMindepNormalGamma_reg_known_vcov()`:

```text
dat (filtered)
   │
   ├─ model_setup(form, dat)              → design  (y, D, group, W, ranks, …)
   │
   ├─ Prior_Setup_GLMM(form, dat, …)      → ps      (pop.* / group.* priors + design)
   │     └─ calls model_setup() internally
   │
   ├─ pfamily_list(ps)                    → pf      (Block~2 dNormal pfamilies)
   │
   └─ dGamma_list(ps, …)                  → prior_list  (Block~1 per-group dGamma)
         │
         └─ rLMMindepNormalGamma_reg_known_vcov(
              y, D, group, W,
              prior_list   = prior_list,
              pfamily_list = pf, …)
```

Ex_13b also calls `model_setup()` once *before* filtering to full-rank schools (`design_all`), then again on the filtered data. The sampler consumes the second `design` plus `pf` / `disp_pf_list`.

## Dedicated docs

| Function | Doc | Role in Ex_13b |
|----------|-----|----------------|
| `model_setup` | [PREFLIGHT_model_setup.md](PREFLIGHT_model_setup.md) | Design matrices + reference `lmer` + rank filter |
| `Prior_Setup_GLMM` | [PREFLIGHT_Prior_Setup_GLMM.md](PREFLIGHT_Prior_Setup_GLMM.md) | Calibrate Block~1/2 priors (calls `model_setup`) |
| `pfamily_list` | [PREFLIGHT_pfamily_list.md](PREFLIGHT_pfamily_list.md) | Wrap Block~2 as `dNormal` pfamilies |
| `dGamma_list` | [PREFLIGHT_dGamma_list.md](PREFLIGHT_dGamma_list.md) | Wrap Block~1 as per-group `dGamma` pfamilies |

Related math (not call trees): [DGAMMA_LIST_MARGINAL_AND_BOUNDS.md](DGAMMA_LIST_MARGINAL_AND_BOUNDS.md), [BLOCK_GIBBS_ERGODICITY_ING.md](BLOCK_GIBBS_ERGODICITY_ING.md), [DESIGN_GLMBAYESCORE_BRIDGE.md](DESIGN_GLMBAYESCORE_BRIDGE.md).

## Migration summary (glmbayesCore)

| Move first | Depends on | Already in glmbayesCore |
|------------|------------|-------------------------|
| `model_setup` + lme4 design / identifiability stack | `lme4`, `Matrix`, `reformulas`; Suggests `glmmTMB` | Nothing on this path |
| `Prior_Setup_GLMM` + ING / pwt helpers | Full `model_setup` tree | `compute_gaussian_prior` (one call) |
| `pfamily_list.Prior_Setup_GLMM` | Prior container fields | `dNormal`, `dIndependent_Normal_Gamma` |
| `dGamma_list.Prior_Setup_GLMM` | Per-group `group.ing_prior` from setup | `dGamma` |

Practical order if migrating: **(1) model_setup**, **(2) Prior_Setup_GLMM**, **(3) thin list methods**. The list methods are small constructors; almost all algebra lives in setup.

## Export status (summary)

Legend used in the four dedicated docs:

| Tag | Meaning |
|-----|---------|
| **Exported** | In `lmebayesCore` `NAMESPACE` today |
| **Internal** | `@noRd` / `@keywords internal` (callable via `:::`) |
| **External** | Exported from another package (`glmbayesCore`, `lme4`, …) |
| **Recommend export** | Should become a documented export if/when this stack moves to `glmbayesCore` (or for advanced lmebayesCore users) |
| **Keep internal** | Plumbing; do not export |

### Already exported on the Ex_13b preflight path

| Symbol | Notes |
|--------|-------|
| `model_setup`, `print.model_setup` | Design + ref fit + ranks |
| `check_identifiability` | Called by `model_setup`; also usable alone |
| `Prior_Setup_GLMM`, `print.Prior_Setup_GLMM` | Prior calibration |
| `pfamily_list`, `pfamily_list.Prior_Setup_GLMM` | Block~2 pfamilies |
| `dGamma_list`, `dGamma_list.Prior_Setup_GLMM` | Block~1 pfamilies |

Related exports **not** on this preflight path but used immediately after: `rLMMindepNormalGamma_reg_known_vcov`, `priors_from_pfamily_list`, `matrix_args_lmm`, `build_mu_all`.

### Strongly recommend exporting (today internal)

These are the main reusable seams if migrating to `glmbayesCore` or supporting design-only workflows:

| Symbol | Why |
|--------|-----|
| `extract_re_hyper_matrices` | Design-only API (`y`/`D`/`W`/`group`) without fitting |
| `extract_mer_variance_components` | Ψ + residual σ² from a merMod (also `extract_lmer_variance_components` alias) |
| `get_lme4_components` | Low-level `lFormula` → X/Z/flist (advanced / debugging) |

### Recommend exporting only with a glmmTMB / calibration package surface

| Symbol | Why |
|--------|-----|
| `.lmebayes_fit_glmmtmb_reference` | Rename without leading `.` if exported; shared by `model_setup` + Prior_Setup |
| `extract_glmmtmb_variance_components` | Parallel to mer VC extract |
| `.lmebayes_glmmtmb_group_sigma2` | Per-group σ² from TMB disp |

### Keep internal (representative)

Formula guards, Z-label maps, weight/offset normalizers, pwt resolvers, ING window helpers, α-target search, Part VI glm inputs — see per-function docs. Exporting them would freeze a large private API for little user benefit.

### Already exported in glmbayesCore (called from this stack)

| Symbol | Called from |
|--------|-------------|
| `compute_gaussian_prior` | Prior_Setup per-group ING |
| `dNormal`, `dIndependent_Normal_Gamma` | `pfamily_list` |
| `dGamma` | `dGamma_list` |

Details and full call-path tables: the four `PREFLIGHT_*.md` files below.
