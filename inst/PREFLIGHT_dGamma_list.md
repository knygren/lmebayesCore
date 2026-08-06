# `dGamma_list` — Block~1 per-group `dGamma` from `Prior_Setup_GLMM`

**Export:** generic `dGamma_list()` + `dGamma_list.Prior_Setup_GLMM()`  
**Source:** `R/dGamma_list.R`, `R/dGamma_list_Prior_Setup_GLMM.R`  
**Ex_13b role:** `disp_pf_list <- dGamma_list(ps, max_disp_perc_measurement = 0.8)` → `prior_list` for the sampler.  
**Index:** [PREFLIGHT_Ex13b_OVERVIEW.md](PREFLIGHT_Ex13b_OVERVIEW.md)  
**Depends on:** [PREFLIGHT_Prior_Setup_GLMM.md](PREFLIGHT_Prior_Setup_GLMM.md)  
**Math:** [DGAMMA_LIST_MARGINAL_AND_BOUNDS.md](DGAMMA_LIST_MARGINAL_AND_BOUNDS.md) (Parts I–VI; Part VI Omega_j is now inside Prior_Setup)

---

## 1. Purpose

Convert precomputed per-group Block~1 measurement-dispersion calibration (`object$group.ing_prior`) into a named list of `glmbayesCore::dGamma()` pfamilies (one per group level).

Like `pfamily_list`, this method is thin. Heavy work (within-group GLM, Part VI Omega_j, `compute_gaussian_prior`, α-target pwt search) already ran in `Prior_Setup_GLMM` when `dispformula = ~group`.

Optional: recompute truncation windows via `max_disp_perc_measurement` without redoing the Gamma shape/rate calibration.

---

## 2. Signature

```r
dGamma_list(object, ...)                                   # generic
dGamma_list.Prior_Setup_GLMM(object, max_disp_perc_measurement = NULL, ...)
```

| Arg | Expectation |
|-----|-------------|
| `object` | `"Prior_Setup_GLMM"`, **gaussian**, with **grouped** `group.ing_prior` (`dispformula = ~group`) |
| `max_disp_perc_measurement` | `NULL` (reuse stored bounds), or scalar / length-`J` in `(0.5, 1)` to rebuild windows |
| `...` | Ignored (legacy args such as old `disp_upper_anchor` are no-ops if passed) |

---

## 3. Call hierarchy (method)

```text
dGamma_list(object, ...)
└─ dGamma_list.Prior_Setup_GLMM
   ├─ stop unless object$family$family == "gaussian"
   ├─ .lmebayes_ing_prior_is_grouped(object$group.ing_prior)
   │    [defined in R/Prior_Setup_GLMM.R; also used by print]
   ├─ if max_disp_perc_measurement != NULL:
   │    └─ .lmebayes_expand_scalar_or_vector(..., range = c(0.5, 1))
   └─ for each group level lev:
        g ← group.ing_prior[[lev]]
        ├─ if mdp differs from g$max_disp_perc or bounds missing:
        │    shape_post = g$shape_ING + g$n_j/2
        │    rate_post  = g$sigma2_hat * (shape_post − 1)
        │    └─ .lmebayes_ing_prior_quantile_window(shape_post, rate_post, mdp)
        │         └─ stats::qgamma
        └─ ★ glmbayesCore::dGamma(
             shape = g$shape_ING,
             rate  = g$rate,
             beta  = 0×1 "(Intercept)",
             Inv_Dispersion = TRUE,
             max_disp_perc, disp_lower, disp_upper)
```

---

## 4. Upstream hierarchy (fills `group.ing_prior`; not called by `dGamma_list`)

This is what Ex_13b’s `Prior_Setup_GLMM(..., dispformula = ~school_id)` runs so that a plain `dGamma_list(ps)` works:

```text
Prior_Setup_GLMM  [per-group branch]
└─ .lmebayes_calibrate_ing_prior_measurement_group()
   ├─ .lmebayes_ing_prior_measurement_group_glm_inputs()
   │   ├─ stats::model.frame / model.matrix / glm.fit / vcov
   │   └─ .lmebayes_block_formula_prior_mu()
   ├─ Sigma_j shrinkage from τ_sd
   ├─ Omega_j ← W_j %*% pop.prior_list[[k]]$Sigma %*% t(W_j)   [Part VI]
   ├─ .lmebayes_compute_ing_prior_cal_from_sigma()
   │   └─ ★ glmbayesCore::compute_gaussian_prior()
   ├─ .ing_stop_if_prior_exceeds_data()
   ├─ .lmebayes_ing_prior_quantile_window()
   └─ .lmebayes_ing_prior_list_from_cal() → group.ing_prior[[lev]]

[+ optional] .lmebayes_calibrate_pwt_measurement_group()
   └─ re-calibrates group.ing_prior after α_target pwt search
```

Ex_13b formerly hand-rolled Part VI via `.lmebayes_ing_prior_measurement_group_glm_inputs` / `.lmebayes_compute_ing_prior_cal_from_sigma`; that block is now commented out — production path is `dGamma_list(ps, …)`.

---

## 5. Return and attributes

Named list of `"pfamily"` (`dGamma`) keyed by group levels, plus:

| Attribute | Contents |
|-----------|----------|
| `"window_diagnostics"` | data.frame: group, n_j, sigma2_hat, shape_ING, rate, max_disp_perc, bounds |
| `"measurement_prior_group"` | Flattened named vectors: `shape_group`, `rate_group`, `disp_lower_group`, `disp_upper_group` |
| `"group.dispersion.fit"` | `object$group.dispersion.fit` (glmmTMB ref when present) |
| `"calibration_source"` | `"lme4"` / `"glmmTMB"` |

Ex_13b reads `attr(disp_pf_list, "measurement_prior_group")` for diagnostics and passes the list itself as `prior_list`.

---

## 6. `Prior_Setup_GLMM` fields used (by the method)

| Field | Role |
|-------|------|
| `family$family` | Must be `"gaussian"` |
| `group.ing_prior` (grouped) | Per-level shape_ING, rate, sigma2_hat, n_j, bounds |
| `design$group` / `group_name` | Level names / errors |
| `group.dispersion.fit` | Echoed as attribute |
| `calibration_source` | Echoed as attribute |

Not read at wrap time: `block_formula`, `data`, `pop.prior_list`, `group.Sigma` (those feed setup earlier).

---

## 7. Helpers

| Symbol | Exclusive? | Purpose |
|--------|------------|---------|
| `.lmebayes_ing_prior_is_grouped` | No | Detect grouped prior shape |
| `.lmebayes_expand_scalar_or_vector` | No | Expand mdp vector |
| `.lmebayes_ing_prior_quantile_window` | No | Rebuild bounds |

No helpers exist solely for this method.

---

## 8. Export status of symbols on this call path

### Already exported (`lmebayesCore` NAMESPACE)

| Symbol | Role |
|--------|------|
| `dGamma_list` | Generic |
| `dGamma_list.Prior_Setup_GLMM` | S3 method |

### External exports used (by this method)

| Symbol | Package |
|--------|---------|
| `glmbayesCore::dGamma` | **Exported** in glmbayesCore |
| `stats::qgamma` | Via `.lmebayes_ing_prior_quantile_window` when rebuilding bounds |

### External exports used **upstream** (Prior_Setup → `group.ing_prior`, not by this method)

| Symbol | Package |
|--------|---------|
| `glmbayesCore::compute_gaussian_prior` | **Exported** in glmbayesCore |

### Internal today — **keep internal**

| Symbol | Recommendation |
|--------|----------------|
| `.lmebayes_ing_prior_is_grouped` | **Keep internal.** Structural check for print + this method. |
| `.lmebayes_expand_scalar_or_vector` | **Keep internal.** Shared vector expand helper. |
| `.lmebayes_ing_prior_quantile_window` | **Keep internal.** Window math; public knob is `max_disp_perc_measurement` on `dGamma_list`. |

### Related exports (not called by this method)

| Symbol | Notes |
|--------|-------|
| `Prior_Setup_GLMM` | Fills `group.ing_prior` |
| `rLMMindepNormalGamma_reg_known_vcov` (and siblings) | **Already exported**; consume the returned list as `prior_list` |

No additional lmebayesCore exports are recommended solely for the `dGamma_list` wrap step — the public surface is already the generic + method.

---

## 9. Migrating to glmbayesCore

### Minimal move (wrap only)

| Piece | Export? | Notes |
|-------|---------|-------|
| Generic + method | **Keep exported** | Thin |
| Grouped `group.ing_prior` contract | — | Must exist somehow |
| Shared window helpers above | **Keep internal** | Also used by Prior_Setup / pwt |

### Already in glmbayesCore (exported there)

- `dGamma`
- `compute_gaussian_prior` (upstream in Prior_Setup, not in this method)

### Realistic migration bundle

Moving `dGamma_list` **without** the per-group calibration path is useless for Ex_13b-style workflows. Migrate together:

1. `Prior_Setup_GLMM` per-group branch ([PREFLIGHT_Prior_Setup_GLMM.md](PREFLIGHT_Prior_Setup_GLMM.md) §7 Layers C/D) — helpers stay internal
2. This method as the **exported** constructor
3. Docs: `DGAMMA_LIST_MARGINAL_AND_BOUNDS.md`

### Attribute consumers (optional follow-ons)

- `rLMMindepNormalGamma_reg_*` (`prior_list`) — already exported
- `lmerb` / `glmerb` reuse of `attr(*, "group.dispersion.fit")`

---

## 10. Ex_13b usage

```r
disp_pf_list <- dGamma_list(ps, max_disp_perc_measurement = 0.8)
prior_list   <- disp_pf_list
mpg <- attr(disp_pf_list, "measurement_prior_group")

fit <- rLMMindepNormalGamma_reg_known_vcov(
  n = 1000L,
  y = design$y, D = design$D, group = grp, W = design$W,
  prior_list   = prior_list,
  pfamily_list = pf,   # from pfamily_list(ps)
  …
)
```
