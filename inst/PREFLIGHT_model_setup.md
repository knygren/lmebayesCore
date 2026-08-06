# `model_setup` — design, reference fit, identifiability

**Export:** `model_setup()` / `print.model_setup()`  
**Source:** `R/model_setup.R` (plus helpers in `R/lme4_design_utilities.R`, `R/check_identifiability.R`, `R/glmmtmb_reference_helpers.R`)  
**Ex_13b role:** Build `design` (`y`, `D`, `group`, `W`, ranks). Also used for the full-rank school filter before the final call.  
**Index:** [PREFLIGHT_Ex13b_OVERVIEW.md](PREFLIGHT_Ex13b_OVERVIEW.md)

---

## 1. Purpose

Turn a single-grouping-factor mixed formula + data into:

1. **Sampler design** — `y`, `D` (within-group loadings), `group`, `W` (hyper-design list), `weights`, `offset`
2. **Reference classical fit** — `lmer` / `glmer` (optional `glmmTMB` when `dispformula = ~group`)
3. **Identifiability diagnostics** — `groupef.rank` / `groupef.estimable` / `popef.rank_ok`

Does **not** calibrate Bayesian priors (that is `Prior_Setup_GLMM`).

---

## 2. Signature (main formals)

```r
model_setup(
  formula, data = NULL, family = gaussian(),
  REML = TRUE, control = NULL, start = NULL, verbose = 0L,
  subset, weights, na.action, offset, contrasts = NULL,
  devFunOnly = FALSE, fit_mer = TRUE, dispformula = ~1, ...
)
```

`subset` / `na.action` / `contrasts` / `weights` / `offset` are shared between design extraction (`lFormula`) and the reference fit so rows and factor coding match.

---

## 3. Call hierarchy

External packages in **bold**.

```text
model_setup()
│
├─ .lmebayes_normalize_family()
├─ lme4::lmerControl()                    [default when gaussian & control NULL]
│
├─ extract_re_hyper_matrices()            [R/lme4_design_utilities.R]
│   ├─ is_single_factor_model()
│   │   ├─ reformulas::findbars()
│   │   └─ lme4::lFormula()
│   ├─ .lmebayes_validate_uncorrelated_re_formula()
│   │   ├─ is_fixed_effects_only()
│   │   │   ├─ reformulas::findbars()
│   │   │   └─ lme4::lFormula()
│   │   └─ lme4::lFormula()               [cnms: reject correlated RE]
│   ├─ lme4::lFormula()                   [main parse → fr, X, reTrms]
│   ├─ classify_lme4_fixed_columns()      [level-1 vs level-2 FE]
│   ├─ classify_crosslevel_re_moderation()[slope means + moderator:slope]
│   ├─ get_lme4_components()
│   │   ├─ lme4::lFormula()
│   │   ├─ lme4::mkLmerDevfun()           [expose pp$Zt]
│   │   ├─ Matrix::t()
│   │   └─ .lme4_label_Z_random_sparse()
│   │       ├─ .lme4_Z_random_rownames() / _colnames()
│   │       └─ .lme4_Z_random_column_map() / _row_map()
│   ├─ extract_lme4_fixed_group_matrix()
│   ├─ extract_re_Z_obs()                 [sparse Z → dense D (n × p_re)]
│   ├─ stats::model.response(fr)
│   └─ .lmebayes_weights_offset_from_frame()
│       ├─ stats::model.weights / model.offset
│       ├─ .lmebayes_normalize_weights()
│       └─ .lmebayes_normalize_offset()
│
├─ stats::terms()
├─ .lmebayes_dispformula_kind()           [R/glmmtmb_reference_helpers.R]
│
├─ [if fit_mer]
│   ├─ .lmebayes_mer_optional_args()
│   ├─ lme4::lmer()  OR  lme4::glmer()
│   ├─ lme4::isSingular() → message()
│   ├─ stats::fitted()                    [n_obs guard vs design$y]
│   ├─ extract_mer_variance_components()
│   │   ├─ lme4::VarCorr()
│   │   └─ lme4::getME(fit, "sigma")^2    [gaussian residual]
│   └─ [gaussian & dispformula = ~group]
│       └─ .lmebayes_fit_glmmtmb_reference()
│           └─ glmmTMB::glmmTMB()
│
└─ check_identifiability()                [export, R/check_identifiability.R]
    ├─ .lmebayes_normalize_family()
    ├─ Matrix::rankMatrix(D_j)            → groupef.rank
    ├─ .lmebayes_block_glm_estimable()
    │   ├─ .lmebayes_glm_estimable_precheck()
    │   └─ stats::glm() / vcov()          [binomial/poisson/Gamma]
    └─ Matrix::rankMatrix(W[[k]][estimable,]) → popef.*
```

### Stages (same tree, by concern)

| Stage | What runs |
|-------|-----------|
| Design / model frame | `extract_re_hyper_matrices` and children |
| Reference fit | `lmer`/`glmer`, VC extract, optional `glmmTMB` |
| Identifiability | `check_identifiability` |

---

## 4. Return (class `"model_setup"`)

| Field | Role |
|-------|------|
| `y`, `weights`, `offset` | Length-`n`, aligned after subset/NA |
| `D` | Within-group design (`n × p_re`) |
| `group`, `group_name`, `groupef.names` | Grouping factor + RE coef names |
| `W` | Named list of hyper-design matrices (one per RE coef) |
| `popef.moderation` | Cross-level moderation metadata |
| `lmer` / `glmer`, `Psi`, `dispersion`, `varcorr` | Reference fit (if `fit_mer`) |
| `glmmTMB_fit` | Additive per-group-dispersion fit, or `NULL` |
| `groupef.rank`, `groupef.estimable`, `groupef.glm_check` | Level-1 |
| `popef.rank`, `popef.deficient`, `popef.rank_ok` | Level-2 |

---

## 5. Export status of symbols on this call path

### Already exported (`lmebayesCore` NAMESPACE)

| Symbol | Role on path |
|--------|----------------|
| `model_setup` | Entry point |
| `print.model_setup` | S3 printer |
| `check_identifiability` | Always called at end of `model_setup` |

### External exports used (not lmebayesCore)

| Symbol | Package |
|--------|---------|
| `lme4::lFormula`, `lmer`, `glmer`, `lmerControl`, `VarCorr`, `getME`, `isSingular`, `mkLmerDevfun` | **lme4** |
| `reformulas::findbars` | **reformulas** |
| `Matrix::rankMatrix`, `Matrix::t` | **Matrix** |
| `glmmTMB::glmmTMB` | **glmmTMB** (Suggests; only if `dispformula = ~group`) |
| `stats::terms`, `model.response`, `model.weights`, `model.offset`, `fitted`, `glm`, `vcov`, … | **stats** |

No **glmbayesCore** exports on this path.

### Internal today — **recommend export** (migration / advanced API)

| Symbol | Current tag | Recommendation |
|--------|-------------|----------------|
| `extract_re_hyper_matrices` | `@keywords internal` | **Recommend export.** Design-only seam (`fit_mer`-free); Ex_13b-style tooling and glmbayesCore migration both need it as a public “formula → y/D/W/group” function. |
| `extract_mer_variance_components` | `@keywords internal` | **Recommend export.** Stable “merMod → Ψ + dispersion” helper; useful without re-running full `model_setup`. |
| `extract_lmer_variance_components` | `@keywords internal` | Thin alias of the mer extractor; export with the same docs or collapse into one name. |
| `get_lme4_components` | `@keywords internal` | **Optional export** for advanced users / debugging (X, sparse Z, maps). Keep `extract_re_hyper_matrices` as the usual public surface. |

### Internal today — **recommend export only if** shipping a glmmTMB surface

| Symbol | Recommendation |
|--------|----------------|
| `.lmebayes_fit_glmmtmb_reference` | Export under a non-dot name (e.g. `fit_glmmtmb_reference`) if glmmTMB calibration is a first-class glmbayesCore feature. |
| `.lmebayes_dispformula_kind` | Keep internal unless documenting a public dispformula grammar helper. |

### Internal today — **keep internal**

| Symbol | Why |
|--------|-----|
| `.lmebayes_normalize_family` | Tiny coerce; not a user API |
| `.lmebayes_mer_optional_args` | Call packing |
| `.lmebayes_normalize_weights` / `_offset` / `_weights_offset_from_frame` | Frame plumbing |
| `is_single_factor_model`, `is_fixed_effects_only` | Formula guards |
| `.lmebayes_validate_uncorrelated_re_formula` | Formula guard |
| `classify_lme4_fixed_columns`, `classify_crosslevel_re_moderation` | Classification helpers |
| `extract_lme4_fixed_group_matrix`, `extract_re_Z_obs` | Steps inside `extract_re_hyper_matrices` |
| `.lme4_Z_random_*`, `.lme4_label_Z_random_sparse` | Sparse Z bookkeeping |
| `.lmebayes_block_glm_estimable`, `.lmebayes_glm_estimable_precheck` | Inside `check_identifiability` |

---

## 6. Co-located symbols **not** called by `model_setup`

These live nearby but are used by `Prior_Setup_GLMM` / samplers:

- `.lmebayes_mer_convergence_issues` — **keep internal**
- `.lmebayes_glmmtmb_convergence_issues`, `.lmebayes_reference_fixef` / `_vcov` — **keep internal** (dispatch helpers)
- `extract_glmmtmb_variance_components`, `.lmebayes_glmmtmb_group_sigma2` — **recommend export with glmmTMB surface** (see §5)
- `.lmebayes_stop_if_nondefault_weights_offset` — **keep internal** (Phase~1 gate)
- `lmerb_default_vcov_formula` — **optional export** if users need the stripped calibration formula outside lmerb

Migrate them with the same modules if you want a complete mixed-design package surface.

---

## 7. Migrating to glmbayesCore

### Already in glmbayesCore?

**Nothing** on this call tree. `glmbayesCore::Prior_Setup` is flat-GLM prior setup (conceptual cousin only).

### Public API to expose after the move

| Export | Notes |
|--------|-------|
| `model_setup`, `print.model_setup` | Already exported here |
| `check_identifiability` | Already exported here |
| `extract_re_hyper_matrices` | **Newly public** (see §5) |
| `extract_mer_variance_components` | **Newly public** (see §5) |
| Optional: `get_lme4_components`, glmmTMB helpers | Advanced / Suggests surface |

### Must move (on-path internals)

All symbols listed in §5 “keep internal” plus the recommend-export candidates — they travel with the modules even if not all become exports.

### External deps that travel with the move

`lme4`, `Matrix`, `reformulas`, `stats`; Suggests / `requireNamespace("glmmTMB")`.

### Suggested packaging split inside glmbayesCore

1. **Design-only** (`fit_mer = FALSE` or bare `extract_re_hyper_matrices`): + `check_identifiability`  
2. **Reference fit** layer: merMod + optional glmmTMB  
3. Keep formula constraints (one group, uncorrelated RE, level-2 FE rules) documented as the GLMM design dialect shared by lmebayes.

---

## 8. Ex_13b usage

```r
design_all <- model_setup(form_lmer, data = dat)          # filter ranks
# … keep full_rank_schools …
design <- model_setup(form_lmer, data = dat)              # final design
# sampler uses design$y, design$D, design$group, design$W
```
