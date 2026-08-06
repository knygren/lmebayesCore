# `Prior_Setup_GLMM` — Block~1 / Block~2 prior calibration

**Export:** `Prior_Setup_GLMM()` / `print.Prior_Setup_GLMM()`  
**Source:** `R/Prior_Setup_GLMM.R` (helpers in `R/mixed_rmerb_helpers.R`, `R/pwt_measurement_group_calibration.R`, `R/glmmtmb_reference_helpers.R`, `R/ing_prior_guard.R`, `R/two_block_ergodicity_ing_marginal.R`)  
**Ex_13b role:** Build `ps` with `dispformula = ~school_id`, `group.alpha_target` pwt search, Part VI `group.ing_prior`.  
**Index:** [PREFLIGHT_Ex13b_OVERVIEW.md](PREFLIGHT_Ex13b_OVERVIEW.md)  
**Depends on:** [PREFLIGHT_model_setup.md](PREFLIGHT_model_setup.md)

---

## 1. Purpose

Calibrate Bayesian priors for the two-block mixed sampler:

| Block | Quantity | Typical Ex_13b product |
|-------|----------|------------------------|
| Block~2 | Population FE priors + τ² | `pop.prior_list`, `pop.ing_prior` → `pfamily_list(ps)` |
| Block~1 | Measurement σ² (pooled or per-group) | `group.ing_prior` → `dGamma_list(ps)` |

Always embeds a full `model_setup` object as `ps$design`.

---

## 2. Signature (main formals)

```r
Prior_Setup_GLMM(
  formula, data, family = gaussian(),
  REML = TRUE, control = NULL, start = NULL, verbose = 0L,
  subset, weights, na.action, offset, contrasts = NULL,
  dispformula = ~1,
  pop.pwt = NULL, pop.pwt_default_low = 0.01, pop.pwt_default_high = 0.05,
  pop.nprior = NULL, pop.sd = NULL, pop.mu = NULL,
  pop.dispersion.pwt = NULL, pop.dispersion.nprior = NULL,
  pop.max_disp_perc = 0.99,
  pop.intercept_source = c("null_model", "full_model"),
  pop.effects_source   = c("null_effects", "full_model"),
  group.dispersion = NULL,
  group.dispersion.pwt = NULL, group.dispersion.nprior = NULL,
  group.alpha_target = 0.01,
  group.max_disp_perc = 0.99
)
```

Ex_13b uses: `pop.pwt = 0.01`, `dispformula = ~school_id`, `group.max_disp_perc = 0.8`, `group.dispersion.pwt = 0.1` (floor), `group.alpha_target = 0.01`.

---

## 3. Call hierarchy

```text
Prior_Setup_GLMM()
│
├─ [arg checks; family normalize; default lmerControl for gaussian]
│
├─ model_setup(...)                         ★ full tree: see PREFLIGHT_model_setup.md
│
├─ .lmebayes_dispformula_kind()
├─ .lmebayes_mer_convergence_issues()
├─ .lmebayes_validate_group_dispersion()
│
├─ [choose calibration_source]
│   ├─ if gaussian & dispformula = ~group:
│   │   ├─ design$glmmTMB_fit (from model_setup) or refit
│   │   ├─ .lmebayes_glmmtmb_convergence_issues()
│   │   ├─ extract_glmmtmb_variance_components()
│   │   └─ .lmebayes_glmmtmb_group_sigma2()
│   └─ else: fit_ref = mer_fit; τ² = design$Psi
│
├─ .lmebayes_reference_fixef() / .lmebayes_reference_vcov()
│
├─ [Block~2 weight / override resolution]
│   ├─ .lmebayes_default_pwt_list()
│   ├─ .lmebayes_resolve_pwt()
│   ├─ .lmebayes_resolve_mu_override() → .lmebayes_resolve_override_list()
│   ├─ .lmebayes_resolve_nprior_override()
│   ├─ .lmebayes_resolve_sd_override()
│   └─ .lmebayes_resolve_disp_prior()       [pop.dispersion.*]
│
├─ [null intercept model] when pop.intercept_source == "null_model"
│   ├─ stats::as.formula(y ~ 1 + (1 | group))
│   ├─ .lmebayes_fit_glmmtmb_reference() OR lme4::lmer/glmer
│   ├─ convergence checks
│   └─ .lmebayes_reference_fixef(null_fit)
│
├─ [inline] build pop.prior_list: mu, Sigma ← V_fe × √((1−w)/w), τ²_k
│
├─ .lmebayes_expand_scalar_or_vector()      [pop.max_disp_perc]
├─ .lmebayes_ing_limiting_posterior_window()→ stats::qgamma
│     [pop.ing_prior shape/rate/bounds; same algebra as compute_gaussian_prior k=1]
│
├─ [Gaussian Block~1 measurement]
│   ├─ .lmebayes_resolve_measurement_disp_prior()          [pooled]
│   ├─ [dispformula = ~1]
│   │   └─ .lmebayes_calibrate_ing_prior_measurement()
│   │       └─ .lmebayes_ing_prior_quantile_window()
│   │
│   ├─ .lmebayes_block_formula_from_re()
│   ├─ .lmebayes_resolve_measurement_disp_prior_group()
│   ├─ .lmebayes_expand_scalar_or_vector()                 [group.max_disp_perc]
│   │
│   ├─ [dispformula = ~group]  ★ Ex_13b path
│   │   └─ .lmebayes_calibrate_ing_prior_measurement_group()
│   │       ├─ .lmebayes_ing_prior_measurement_group_glm_inputs()
│   │       │   ├─ stats::model.frame / model.matrix / model.response
│   │       │   ├─ stats::glm.fit / vcov
│   │       │   └─ .lmebayes_block_formula_prior_mu() → stats::lm
│   │       ├─ Omega_j ← W_j %*% pop.Sigma_k %*% t(W_j)   [Part VI]
│   │       ├─ .lmebayes_compute_ing_prior_cal_from_sigma()
│   │       │   └─ ★ glmbayesCore::compute_gaussian_prior()
│   │       ├─ .ing_stop_if_prior_exceeds_data()
│   │       │   └─ .ing_n_prior_from_shape()
│   │       ├─ .lmebayes_ing_prior_quantile_window()
│   │       └─ .lmebayes_ing_prior_list_from_cal()
│   │
│   ├─ .lmebayes_print_ing_prior_measurement_group_compare()  [dev table]
│   │
│   └─ [group.alpha_target]  ★ Ex_13b path
│       ├─ .lmebayes_calibrate_pwt_measurement_group()
│       │   ├─ .lmebayes_reference_fixef()
│       │   ├─ stats::rnorm / uniroot
│       │   ├─ .lmebayes_ing_prior_quantile_window()
│       │   └─ .two_block_truncated_omega_moments()
│       │         [R/two_block_ergodicity_ing_marginal.R]
│       ├─ re-resolve meas pwt + re-calibrate group.ing_prior
│       └─ store group.pwt_calibration
│
└─ [optional] .lmebayes_pin_ing_prior_measurement_group()
    when user overrides group.dispersion
```

### Only direct `glmbayesCore` runtime call

| Call | Where | When |
|------|--------|------|
| `glmbayesCore::compute_gaussian_prior` | `.lmebayes_compute_ing_prior_cal_from_sigma` | Gaussian + per-group `dispformula` |

Block~2 `pop.ing_prior` **reimplements** the same shape/rate algebra locally (no `compute_gaussian_prior` call).

---

## 4. Return (class `"Prior_Setup_GLMM"`)

| Field | Role |
|-------|------|
| `design` | Full `"model_setup"` object |
| `mer_fit`, `fit_ref`, `calibration_source` | Classical references |
| `pop.pwt`, `pop.nprior`, `pop.sd`, `pop.mu` | Resolved Block~2 strength / mean |
| `pop.prior_list` | Per-RE `{mu, Sigma, dispersion=τ²}` |
| `pop.ing_prior` | Per-RE τ² Gamma window |
| `pop.dispersion.pwt` / `pop.dispersion.nprior` | τ² prior strength |
| `group.dispersion`, `group.Sigma`, `group.tau_sd` | Measurement / RE variances |
| `group.dispersion.pwt` / `group.dispersion.nprior` | Block~1 strength |
| `group.ing_prior` | Pooled list **or** named per-group lists (Ex_13b) |
| `group.alpha_target`, `group.pwt_calibration` | Ellipsoid pwt search |
| `group.max_disp_perc`, `block_formula` | Window / within-group formula |

---

## 5. Reuse vs Prior_Setup-only

| Reused from `model_setup` | Prior_Setup-specific |
|---------------------------|----------------------|
| Design `y/D/W/group` | `pop.pwt` / `pop.nprior` / `pop.sd` duality |
| merMod + optional `glmmTMB_fit` | Null intercept model for prior mean |
| `Psi`, pooled `dispersion` | Block~2 `Sigma` scaling from `vcov` |
| Rank diagnostics (reported) | `pop.ing_prior`, `group.ing_prior` |
| | Part VI Omega_j + `compute_gaussian_prior` |
| | `group.alpha_target` Monte Carlo / `uniroot` |

---

## 6. Export status of symbols on this call path

### Already exported (`lmebayesCore` NAMESPACE)

| Symbol | Role on path |
|--------|----------------|
| `Prior_Setup_GLMM` | Entry point |
| `print.Prior_Setup_GLMM` | S3 printer |
| `model_setup` (+ `print.model_setup`) | Called internally for design / ref fits |
| `check_identifiability` | Via `model_setup` |

Downstream list constructors (next Ex_13b steps, not called *by* Prior_Setup): `pfamily_list`, `dGamma_list` — already exported.

### External exports used

| Symbol | Package | Role |
|--------|---------|------|
| **`glmbayesCore::compute_gaussian_prior`** | glmbayesCore | Per-group Block~1 ING calibration (only glmbayesCore runtime call) |
| `lme4::lmer` / `glmer` / `lmerControl` / … | lme4 | Via `model_setup` + null intercept fit |
| `glmmTMB::glmmTMB`, `VarCorr`, `fixef` | glmmTMB | Per-group dispformula path |
| `stats::qgamma`, `uniroot`, `rnorm`, `glm.fit`, `lm`, … | stats | Windows, α-target search, within-group fits |

### Internal today — **recommend export** (with Prior_Setup / calibration surface)

| Symbol | Recommendation |
|--------|----------------|
| `extract_re_hyper_matrices`, `extract_mer_variance_components` | Same as [PREFLIGHT_model_setup.md](PREFLIGHT_model_setup.md) §5 — Prior_Setup depends on them via `model_setup`. |
| `.lmebayes_fit_glmmtmb_reference` | Export (rename without `.`) if glmmTMB is a public calibration backend. |
| `extract_glmmtmb_variance_components` | Export alongside mer VC extractor for a symmetric ref-fit API. |
| `.lmebayes_glmmtmb_group_sigma2` | Export if users need per-group σ² diagnostics without full Prior_Setup. |

### Internal today — **keep internal** (do not export)

**Prior_Setup-local resolvers** (`R/Prior_Setup_GLMM.R`):  
`.lmebayes_validate_group_dispersion`, `.lmebayes_default_pwt_list`, `.lmebayes_resolve_pwt`, `.lmebayes_resolve_override_list`, `.lmebayes_resolve_mu_override` / `_sd_override` / `_nprior_override`, `.lmebayes_resolve_disp_prior`, `.lmebayes_ing_prior_is_grouped`.

**Measurement / ING helpers** (`R/mixed_rmerb_helpers.R`):  
`.lmebayes_expand_scalar_or_vector`, `.lmebayes_resolve_measurement_disp_prior` / `_group`, `.lmebayes_block_formula_from_re` / `_prior_mu`, `.lmebayes_ing_prior_measurement_group_glm_inputs`, `.lmebayes_compute_ing_prior_cal_from_sigma`, `.lmebayes_ing_prior_list_from_cal`, `.lmebayes_calibrate_ing_prior_measurement` / `_group`, `.lmebayes_ing_prior_quantile_window`, `.lmebayes_ing_limiting_posterior_window`, `.lmebayes_pin_ing_prior_measurement_group*`, `.lmebayes_print_ing_prior_measurement_group_compare`.

**α-target / guards:**  
`.lmebayes_calibrate_pwt_measurement_group`, `.two_block_truncated_omega_moments`, `.ing_stop_if_prior_exceeds_data`, `.ing_n_prior_from_shape`.

**Dispatch / conv:**  
`.lmebayes_mer_convergence_issues`, `.lmebayes_glmmtmb_convergence_issues`, `.lmebayes_reference_fixef` / `_vcov`, `.lmebayes_dispformula_kind`.

Rationale: these encode calibration policy and numerical details; the stable user surface is `Prior_Setup_GLMM` itself (plus optional design/VC extractors above).

### Optional export (only if supporting hand-rolled Part VI demos)

Ex_13b historically called `.lmebayes_ing_prior_measurement_group_glm_inputs` and `.lmebayes_compute_ing_prior_cal_from_sigma` from demo code. Prefer keeping them internal and documenting `dGamma_list(ps)` as the public path (current Ex_13b production path).

---

## 7. Migrating to glmbayesCore

This is **not** a single-file move. Suggested layers (export policy in §6):

### Layer A — already covered by migrating `model_setup`

Entire [PREFLIGHT_model_setup.md](PREFLIGHT_model_setup.md) on-path set, plus often:

- `.lmebayes_mer_convergence_issues` (**keep internal**)
- Remaining `glmmtmb_reference_helpers.R` (export VC / group-σ² helpers if shipping TMB surface)

### Layer B — Prior_Setup-only (`R/Prior_Setup_GLMM.R`)

| Symbol | Export? | Purpose |
|--------|---------|---------|
| `Prior_Setup_GLMM`, `print.Prior_Setup_GLMM` | **Already / keep exported** | Entry + printer |
| Resolvers / validators listed in §6 | **Keep internal** | Policy plumbing |
| Inline null-model + `pop.prior_list` / `pop.ing_prior` builders | **Keep internal** | Core Block~2 body |

### Layer C — measurement / ING (`R/mixed_rmerb_helpers.R` path)

Move with Prior_Setup; **keep internal** (see §6). Public surface remains `Prior_Setup_GLMM` + `dGamma_list`.

### Layer D — pwt α-target + guards

Move with Prior_Setup; **keep internal**.

### Already in glmbayesCore (exported there)

- `compute_gaussian_prior` (used once on the per-group path)
- Conceptual twin: flat `Prior_Setup` (different API; do not confuse)

### Do not need for Prior_Setup alone

Sampler glue in `mixed_rmerb_helpers.R` after the calibration helpers (`priors_from_pfamily_list` — **already exported** but sampler-facing; engine runners, ICM printers, etc.).

---

## 8. Ex_13b usage

```r
ps <- Prior_Setup_GLMM(
  form_lmer, data = dat,
  pop.pwt = 0.01,
  dispformula = ~school_id,
  group.max_disp_perc = 0.8,
  group.dispersion.pwt = 0.1,
  group.alpha_target = 0.01
)
# → pf <- pfamily_list(ps)
# → disp_pf_list <- dGamma_list(ps, max_disp_perc_measurement = 0.8)
```
