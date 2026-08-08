# `pfamily_list` — Block~2 pfamilies from `Prior_Setup_GLMM`

**Export:** generic `pfamily_list()` + `pfamily_list.Prior_Setup_GLMM()`  
**Source:** `R/pfamily_list.R`, `R/pfamily_list_Prior_Setup_GLMM.R`  
**Ex_13b role:** `pf <- pfamily_list(ps)` → Block~2 `dNormal` (known τ²) for every RE component.  
**Index:** [PREFLIGHT_Ex13b_OVERVIEW.md](PREFLIGHT_Ex13b_OVERVIEW.md)  
**Depends on:** [PREFLIGHT_Prior_Setup_GLMM.md](PREFLIGHT_Prior_Setup_GLMM.md)

---

## 1. Purpose

Convert calibrated Block~2 fields on a `"Prior_Setup_GLMM"` object into a named list of `glmbayesCore` **pfamily** objects (one per random-effect coefficient).

This is a thin constructor: almost all algebra already ran inside `Prior_Setup_GLMM`. Downstream samplers (`rLMMindepNormalGamma_reg_*`, `two_block_rNormal_reg`, `lmerb`) take `pfamily_list` as Block~2 prior input.

---

## 2. Signature

```r
pfamily_list(object, ...)                          # generic
pfamily_list.Prior_Setup_GLMM(object, ptypes = "dNormal", ...)
```

| Arg | Expectation |
|-----|-------------|
| `object` | Class `"Prior_Setup_GLMM"` |
| `ptypes` | `"dNormal"` or `"dIndependent_Normal_Gamma"`: scalar (recycled), length-`p_re` character/list, optionally named by RE coef |

Ex_13b uses the default `ptypes = "dNormal"` (known RE variances at REML Ψ).

---

## 3. Call hierarchy

```text
pfamily_list(object, ...)
└─ pfamily_list.Prior_Setup_GLMM
   ├─ validate / expand ptypes → names(object$pop.prior_list)
   ├─ J ← nlevels(object$design$group)
   ├─ local n_prior_for(k)
   │    └─ object$pop.dispersion.nprior[[k]]
   │       else mean(pop.pwt[[k]]) → (w/(1−w))*J
   └─ for each RE name k:
        pl ← object$pop.prior_list[[k]]   # mu, Sigma, dispersion=τ²_k
        ├─ "dNormal"
        │    └─ ★ glmbayesCore::dNormal(mu, Sigma, dispersion)
        └─ "dIndependent_Normal_Gamma"
             ├─ stop if n_prior_k > J
             ├─ prefer object$pop.ing_prior[[k]]
             │    (shape, rate, disp_lower, disp_upper)
             ├─ fallback if NULL:
             │    shape = (n0+1)/2 + p_k/2
             │    rate  = d_k*(n0+p_k−1)/2
             │    └─ .lmebayes_ing_limiting_posterior_window(d_k, J)
             │         └─ stats::qgamma
             └─ ★ glmbayesCore::dIndependent_Normal_Gamma(
                  mu, Sigma, shape, rate, disp_lower, disp_upper)
```

No `compute_gaussian_prior`, no `dGamma`, no design rebuild.

---

## 4. Return

Named `list` of `"pfamily"` objects, keys = `names(object$pop.prior_list)`.  
No list-level attributes.

---

## 5. `Prior_Setup_GLMM` fields used

| Field | Role |
|-------|------|
| `pop.prior_list[[k]]$mu`, `$Sigma`, `$dispersion` | Always |
| `design$group` | `J` for ING balance guard |
| `pop.dispersion.nprior` | Preferred `n0` for ING |
| `pop.pwt` | Fallback `n0` |
| `pop.ing_prior[[k]]` | Preferred ING shape/rate/bounds |

---

## 6. Helpers

| Symbol | Exclusive to this method? | Purpose |
|--------|---------------------------|---------|
| Local `n_prior_for()` | Yes (closure) | Resolve n0 |
| `.lmebayes_ing_limiting_posterior_window` | No (also Prior_Setup) | Fallback τ² window |

---

## 7. Export status of symbols on this call path

### Already exported (`lmebayesCore` NAMESPACE)

| Symbol | Role |
|--------|------|
| `pfamily_list` | Generic |
| `pfamily_list.Prior_Setup_GLMM` | S3 method |

### External exports used

| Symbol | Package |
|--------|---------|
| `glmbayesCore::dNormal` | **Exported** in glmbayesCore |
| `glmbayesCore::dIndependent_Normal_Gamma` | **Exported** in glmbayesCore |
| `stats::qgamma` | Via fallback window helper |

### Internal today — **keep internal**

| Symbol | Recommendation |
|--------|----------------|
| `.lmebayes_ing_limiting_posterior_window` | **Keep internal.** Fallback only when `pop.ing_prior` is missing; normal Ex_13b path uses Prior_Setup-filled windows. Do not export a second public window API. |
| Local `n_prior_for` | Closure; not a package symbol |

### Related exports (not called by this method)

| Symbol | Notes |
|--------|-------|
| `Prior_Setup_GLMM` | Produces `object` |
| `priors_from_pfamily_list` | **Already exported**; consumes a `pfamily_list` result for sampler prior packing — migrate with engines, not required to *build* the list |

---

## 8. Migrating to glmbayesCore

### Minimal move

| Piece | Export? | Notes |
|-------|---------|-------|
| Generic `pfamily_list` + method | **Keep exported** | Tiny |
| Object field contract above | — | Or a thinner “Block~2 prior container” |
| `.lmebayes_ing_limiting_posterior_window` | **Keep internal** | Only needed for ING fallback |

### Already in glmbayesCore (exported there)

- `dNormal`
- `dIndependent_Normal_Gamma`

### Usually migrate **with** Prior_Setup, not alone

Without `Prior_Setup_GLMM` (or an equivalent that fills `pop.prior_list` / `pop.ing_prior`), the method has nothing to wrap. Prefer migrating setup first, then this method as a one-file follow-on.

### Downstream consumers (optional, not required to *build* the list)

- `priors_from_pfamily_list` (**already exported**)
- Two-block / `rLMM*` validators that check pfamily types

---

## 9. Ex_13b usage

```r
pf <- pfamily_list(ps)   # default dNormal for each RE component
# …
fit <- rLMMindepNormalGamma_reg_known_vcov(
  …,
  pfamily_list   = pf,
  dispprior_list = disp_pf_list,   # from dGamma_list
  …
)
```

Related docs: [ADDING_PFAMILY.md](ADDING_PFAMILY.md), [TAU2_ING_FORMULAS.md](TAU2_ING_FORMULAS.md).
