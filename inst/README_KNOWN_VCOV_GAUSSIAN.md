# Fixed dispersion + known variance components: the `known_vcov` Gaussian route

Function-by-function reference for the **exact multivariate-normal** sampling
path used when a linear mixed model has:

- a **fixed** (known, not sampled) observation-level dispersion
  \eqn{\sigma^2} -- scalar (`"fixed"`) or per-group (`"fixed_vector"`), and
- **known variance components** -- every Block~2 hyperparameter uses a
  `dNormal()` pfamily (a fixed \eqn{\tau^2_k}, not sampled from a `dGamma`/ING
  prior).

Under these two conditions the joint posterior of `(fixef, ranef)` is exactly
multivariate normal, so this route never needs an accept/reject envelope,
ICM *iteration*, or a pilot stage -- it computes the exact joint posterior
mean in closed form and calibrates the number of inner Gibbs sweeps from
Theorem~3 eigenvalues computed on that same closed-form structure.

Since the joint posterior is exactly Gaussian, it can also be sampled from
**directly** -- one Cholesky factorization plus `n` independent draws, no
Gibbs sweeps at all. A `sim_method` argument (\S4.4a) selects between that
exact-iid engine (the new default) and the two-block Gibbs engine described
above (\S4.5-\S4.7); see \S4.4a and \S9 for details.

Companion to [R_EXPORT_REACHABILITY.md](R_EXPORT_REACHABILITY.md) (all
exports) and [LMBBLOCK_LMEBAYESCORE_DEPENDENCIES.md](LMBBLOCK_LMEBAYESCORE_DEPENDENCIES.md)
(call-graph style used here).

**Last reviewed:** 2026-07-20, against `lmebayesCore` `R/rLMM_reg.R`,
`R/two_block_measurement_prior.R`, `R/lmebayes_posterior_icm.R`,
`R/two_block_ergodicity.R`, `R/two_block_rNormal_reg.R`,
`R/rLMMNormal_joint_iid.R`, `R/mixed_rmerb_helpers.R`,
`R/lmebayes_reg_route_table.R`, `R/rlmerb.R`/`R/rglmerb.R`, and `src/*.cpp`.

---

## 1. Scope: which `dispersion_ranef` modes land here

`dispersion_ranef` (`lmerb()`/`glmerb()`) resolves to one of five internal
`disp_mode` values via
[`.lmebayes_resolve_dispersion_ranef()`](../R/mixed_rmerb_helpers.R):

| `disp_mode` | Meaning | Sampled? | Route (`REG_ROUTE_TABLE`) |
|---|---|---|---|
| `"none"` | no observation-level dispersion (non-Gaussian family without one) | -- | `glmm_*` |
| **`"fixed"`** | single positive scalar \eqn{\sigma^2}, same for all observations | No | `lmm_fixed_*` -- **this doc** |
| **`"fixed_vector"`** | named length-\eqn{J} numeric vector, one fixed \eqn{\sigma^2_j} per group | No | `lmm_fixed_*` -- **this doc** |
| `"gamma"` | single `dGamma(Inv_Dispersion = TRUE)` pooled prior | Yes | `lmm_gamma_*` |
| `"gamma_list"` | named list of per-group `dGamma()` priors | Yes | `lmm_gamma_*` |

`"fixed"` and `"fixed_vector"` are the two modes with **no** dispersion
sampling; both are routed by
[`.lmebayes_reg_route_key()`](../R/lmebayes_reg_route_table.R) to the same
`lmm_fixed_known` / `lmm_fixed_estimated` table entries (whichever Block~2
`pfamily_list` is *also* all-`dNormal` or not -- see \S2). This doc covers
the `lmm_fixed_known` cell: fixed dispersion **and** all-`dNormal` Block~2,
i.e. `rLMMNormal_reg_known_vcov()`.

`"gamma"`/`"gamma_list"` (sampled dispersion) and the estimated-\eqn{\tau^2}
(`lmm_fixed_estimated`, `rLMMNormal_reg_estimated_vcov()`, pilot-stage) route
are **out of scope** here -- they use ICM *iteration* / conservative rate
plug-ins instead of the exact machinery below.

---

## 2. Headline finding

For fixed dispersion + all-`dNormal` Block~2 there is **no iterative mode
search and no conservative rate bound** anywhere in the call graph:

- The joint posterior mean is obtained by **one linear solve** (dimension =
  total Block~2 hyperparameter count, independent of the number of groups
  \eqn{J}) via Schur-complement elimination of Block~1 -- see
  [`lmerb_posterior_mean()`](#42-lmerb_posterior_mean-exact-closed-form).
- The Theorem~3 convergence rate is **exact** (`convergence_info$method ==
  "exact"`), not a `disp_lower`/`disp_upper` plug-in bound.
- Both computations key off the *same* per-group \eqn{\sigma^2_j} vector
  (scalar broadcast for `"fixed"`, per-group for `"fixed_vector"`), and both
  were fixed in this feature's implementation to correctly index that vector
  instead of assuming a scalar (see `NEWS.md`, "latent bug" entries).
- Because the target is *exactly* Gaussian, `sim_method = "DEFAULT"` skips
  the Gibbs sweep loop entirely: the same closed-form precision matrices
  (`M`, `post_P_j`) that give the exact posterior *mean* are Cholesky-factored
  *once* and every stored draw is one triangular solve away -- see \S4.4a.

---

## 3. Call graph

```
lmerb(dispersion_ranef = <scalar|named vector>, ...)          [lmebayes/R/lmerb.R]
glmerb(family = gaussian(), dispersion_ranef = ..., ...)      [lmebayes/R/glmerb.R]
  |
  v
rlmerb(n, design, prior, dispersion_ranef, ...)                [lmebayesCore/R/rlmerb.R]
rglmerb(..., family = gaussian())                              [lmebayesCore/R/rglmerb.R]
  |
  +-- .lmebayes_resolve_dispersion_ranef()                     [mixed_rmerb_helpers.R]
  |     +-- .lmebayes_resolve_dispersion_ranef_fixed_vector()  [named-vector branch only]
  |
  +-- .lmebayes_block1_prior_list()                            [Block~1 dispersion glue; no longer
  |                                                               builds P -- routed exports derive it]
  |
  v
.lmebayes_run_lmm_engine()                                     [mixed_rmerb_helpers.R]
  +-- .lmebayes_reg_route_key() -> "lmm_fixed_known"            [lmebayes_reg_route_table.R]
  +-- .lmebayes_reg_route_fn()  -> rLMMNormal_reg_known_vcov     [REG_ROUTE_TABLE lookup]
  +-- .lmebayes_matrix_args_lmm()                               [assembles prior_list = list(dispersion = ...);
  |                                                               args$sim_method <- sim_method on this route only]
  |
  v
rLMMNormal_reg_known_vcov(n, y, x, group, x_hyper, prior_list, pfamily_list, sim_method, ...)
  [rLMM_reg.R]  <-- also reachable via rLMMNormal_reg() dispatcher when it
                    detects an all-dNormal pfamily_list
  +-- .rLMM_validate_sim_method()                                ["DEFAULT" | "TWO_BLOCK_GIBBS"]
  |
  +-- sim_method == "DEFAULT"?  ---------------------+  sim_method == "TWO_BLOCK_GIBBS"?
  |                                                   |
  v                                                   v
rLMMNormal_reg_known_vcov_iid()                    rLMMNormal_reg_known_vcov_two_bg()
  [rLMM_reg.R -- \S4.4a]                              [rLMM_reg.R -- \S4.4/4.5, identical body to the
  +-- .rLMM_validate_matrix_inputs(),                             pre-sim_method rLMMNormal_reg_known_vcov()]
  |     .rLMM_validate_fixed_dispersion_prior_list(), +-- .rLMM_validate_matrix_inputs(),
  |     .two_block_validate_pfamily_list()            |     .rLMM_validate_fixed_dispersion_prior_list()
  |     [no P here -- rLMMNormal_joint_iid() below     |     +-- .rLMM_validate_fixed_dispersion_vector()
  |      derives+validates its own P internally]       +-- .two_block_validate_pfamily_list(),
  +-- .rLMMNormal_reg_run_iid()                        |     .two_block_summarize_pfamily_list(),
  |     [prior_list_block1 = list(dispersion, ddef);   |     .rLMM_P_from_pfamily_list()
  |      no P]                                         |     (stop()s here if any component is not
  |     +-- rLMMNormal_joint_iid()  *** exact iid ***  |     dNormal -- wrong route)
  |           [rLMMNormal_joint_iid.R]                 |
  |           re-validates matrix inputs, then          v
  |           .rLMM_P_from_pfamily_list() derives     .rLMMNormal_reg_run()   [shared run body, no pilot]
  |           + injects P (stop()s if caller           +-- .two_block_validate_block1_prior()
  |           supplied P/Sigma)                        +-- .rLMM_icm_at_start()  [EXACT joint posterior mean]
  |     +-- .rLMM_format_v2_out()                      |     +-- .two_block_icm_at_start()
  |                                                      |           +-- .two_block_measurement_prior_list()
  |                                                      |           +-- lmerb_posterior_mean()  *** exact ***
  |                                                      |                 +-- build_mu_all()
  |                                                      +-- .rLMM_calibrate_m_convergence()  [EXACT Theorem~3]
  |                                                      |     +-- two_block_rate_from_pfamily_list() -> ...
  |                                                      |           -> .two_block_rate_inputs(),
  |                                                      |              .two_block_S_P11(), .two_block_gen_eigen()
  |                                                      |     +-- two_block_l_for_tv(method = "theorem3")
  |                                                      |     +-- .two_block_cap_inner_sweeps()
  |                                                      +-- two_block_rNormal_reg(..., m_convergence = ...)
  |                                                      |     [called with a P-free prior_list_block1;
  |                                                      |      two_block_rNormal_reg() derives its own P]
  |                                                      |     +-- .two_block_rNormal_reg_cpp()  [C++ bridge]
  |                                                      +-- .rLMM_format_v2_out()
  |                                                            +-- .two_block_as_staged_names()
  |                                                                                                        |
  +-------------------------- both branches set staged$sim_method_used and return -----------------------+

.lmebayes_run_lmm_engine() (back in rlmerb()/rglmerb())
  +-- .lmebayes_attach_sigma2()   -- "fixed_vector" returns the same named
        length-J vector as both sigma2 and sigma2.mean (never sampled)

rlmerb()/rglmerb() (back in lmebayes' lmerb()/glmerb())
  +-- .lmebayes_add_fixef_summaries()                            [posterior summaries on output]
```

The `rLMMNormal_reg_known_vcov_iid()` leg never touches
`two_block_rNormal_reg()`/`.two_block_rNormal_reg_cpp()` or any compiled
Gibbs-sweep code at all (\S6 describes the `_two_bg` leg's C++ chain only);
`.lmerb_posterior_normal_system()`/`.lmerb_posterior_system_cholesky()`
(\S4.4a) are pure R, reusing the same exact closed-form matrices as
`lmerb_posterior_mean()`.

---

## 4. R-level functions, in call order

### 4.1 Entry points

| Function | File | Role |
|---|---|---|
| `lmerb()` | `lmebayes/R/lmerb.R` | User-facing formula interface (Gaussian only). Narrows the `glmmTMB` reference-fit trigger to `dispersion_mode == "gamma_list"`, so `"fixed"`/`"fixed_vector"` never pull in `glmmTMB`. |
| `glmerb()` | `lmebayes/R/glmerb.R` | User-facing formula interface; `dispersion_ranef` (incl. `"fixed_vector"`) only applies when `family = gaussian()`. |
| [`rlmerb()`](../R/rlmerb.R) | `lmebayesCore/R/rlmerb.R` | Matrix-level entry: `design` (`model_setup`) + `prior` + `dispersion_ranef` + `sim_method`. Resolves dispersion, validates `sim_method` (`.rLMM_validate_sim_method()`), builds `block1_prior` for reporting, calls `.lmebayes_run_lmm_engine()`, attaches ICM table / summaries. Output gains `sim_method_used`. |
| `rglmerb()` | `lmebayesCore/R/rglmerb.R` | Same shape as `rlmerb()` but for `glmerb()`, incl. `sim_method`; for `family = gaussian()` it funnels through the identical `.lmebayes_run_lmm_engine()` route (ignored, i.e. always two-block Gibbs, for non-Gaussian families). |

### 4.2 Dispersion resolution

| Function | File | Role |
|---|---|---|
| `.lmebayes_resolve_dispersion_ranef()` | `mixed_rmerb_helpers.R` | Dispatches `dispersion_ranef` by shape: `NULL`/scalar -> `"fixed"`; named numeric vector (`length > 1`) -> `"fixed_vector"`; `dGamma()` pfamily -> `"gamma"`; named list of pfamilies -> `"gamma_list"`. |
| `.lmebayes_resolve_dispersion_ranef_fixed_vector()` | `mixed_rmerb_helpers.R` | Validates the named vector has length `J`, names exactly matching `levels(design$groups)`, all-positive/finite; reorders to `group_levels`; returns `list(mode = "fixed_vector", dispersion_fix = <named vector>, ...)`. |
| `.lmebayes_block1_prior_list()` | `mixed_rmerb_helpers.R` | Builds `list(dispersion = <scalar-or-vector>, ddef = FALSE)` -- the `prior_list_block1`/`prior_list` shape consumed throughout \S4.3-4.5. No longer builds `P`: every routed export (`rLMMNormal_reg*`, `rGLMM_reg*`, `two_block_rNormal_reg()`, `rLMMNormal_joint_iid()`) derives its own Block~1 precision from `pfamily_list` and rejects a caller-supplied `P`/`Sigma`. |

### 4.3 Routing

| Function | File | Role |
|---|---|---|
| `.lmebayes_reg_route_key()` | `lmebayes_reg_route_table.R` | `disp_mode %in% c("gamma","gamma_list")` -> `lmm_gamma_*`; everything else (incl. `"fixed_vector"`) + Gaussian -> `lmm_fixed_*`; suffix `known`/`estimated` from `pfamily_list`'s dNormal-ness. |
| `REG_ROUTE_TABLE` / `.lmebayes_reg_route_fn()` | `lmebayes_reg_route_table.R` | Declarative map `lmm_fixed_known -> rLMMNormal_reg_known_vcov` (`needs_pilot = FALSE`). |
| `.lmebayes_matrix_args_lmm()` | `mixed_rmerb_helpers.R` | Assembles the flat argument list for the resolved export; for `"fixed"`/`"fixed_vector"` sets `args$prior_list <- list(dispersion = disp_info$dispersion_fix)` (scalar or named vector, passed through unchanged) **and** `args$sim_method <- sim_method` -- only on this branch, since `rLMMNormal_reg_known_vcov()` is the only export with a real `sim_method` dispatch; every other route ignores it. |
| `.lmebayes_run_lmm_engine()` | `mixed_rmerb_helpers.R` | `route_key -> route_fn -> do.call(route_fn, args) -> .lmebayes_attach_sigma2()`. If the callee didn't set `out$sim_method_used` (every route besides `lmm_fixed_known`, which only has a Gibbs engine), fills in `"TWO_BLOCK_GIBBS"`. |
| `.lmebayes_attach_sigma2()` | `mixed_rmerb_helpers.R` | `"fixed_vector"`: returns the same constant named `J`-vector as both `out$sigma2` and `out$sigma2.mean` (nothing sampled, nothing to average). |

### 4.4 `rLMMNormal_reg_known_vcov()` and its validators

| Function | File | Role |
|---|---|---|
| `rLMMNormal_reg()` | `rLMM_reg.R` | Dispatcher: validates inputs once, then re-dispatches (via `match.call()`) to `rLMMNormal_reg_known_vcov()` or `_estimated_vcov()` depending on whether `pf_summary$all_dNormal`. Has a `sim_method` formal that it forwards unchanged. |
| **`rLMMNormal_reg_known_vcov()`** | `rLMM_reg.R` | The route's named export -- now a **thin `sim_method` dispatcher** (via `match.call()` re-dispatch, `.rLMM_validate_sim_method()` first): `"DEFAULT"` -> `rLMMNormal_reg_known_vcov_iid()` (\S4.4a); `"TWO_BLOCK_GIBBS"` -> `rLMMNormal_reg_known_vcov_two_bg()` (below). Neither branch re-validates twice; each callee runs its own full validation. |
| `rLMMNormal_reg_known_vcov_two_bg()` | `rLMM_reg.R` | The pre-`sim_method` `rLMMNormal_reg_known_vcov()` body, factored out unchanged: validates matrix inputs/dispersion/`pfamily_list`, derives `P` from `pfamily_list` (`.rLMM_P_from_pfamily_list()`), **requires all-`dNormal`** (else `stop()`s, pointing at `_estimated_vcov()`/`rLMMNormal_reg()`), then delegates to `.rLMMNormal_reg_run()` (\S4.5-4.7, two-block Gibbs). |
| `.rLMM_validate_matrix_inputs()` | `rLMM_reg.R` | Shape/type checks on `n, y, x, x_hyper, tv_tol, group_name, group`. `re_names`/`group_levels` are no longer separate arguments -- they are always `colnames(x)` (must be unique, non-empty) and `levels(group)` (`group` must be a factor); `group_name` is resolved by the caller via `.lmebayes_resolve_group_name()` (attribute-on-`group` first, then `substitute(group)`) before this function is called, and this function only sanity-checks the resolved value. |
| `.rLMM_validate_fixed_dispersion_prior_list()` | `rLMM_reg.R` | Requires `prior_list$dispersion` and no unexpected fields; delegates numeric validation to the shared vector validator below. |
| `.rLMM_validate_fixed_dispersion_vector()` | `rLMM_reg.R` | Shared scalar-or-length-`J` validator (also used by `rLMMNormal_reg_estimated_vcov()`): accepts a single positive scalar (broadcast) **or** a length-`J` vector, optionally named -- if named, requires an exact set-match to `group_levels` and reorders accordingly. |
| `.two_block_validate_pfamily_list()` | `two_block_rNormal_reg.R` | Structural validation of the Block~2 `pfamily_list` (one pfamily per RE component, dimensions vs. `x_hyper`/`re_names`). |
| `.two_block_summarize_pfamily_list()` | `two_block_rNormal_reg.R` | Computes `pf_summary$all_dNormal` / `any_non_normal` / `ptypes` used to gate `known_vcov` vs. `estimated_vcov` and to tag `convergence_info$method`. |
| `.rLMM_P_from_pfamily_list()` | `two_block_tau2_ref.R` | Derives the Block~2 random-effect prior precision `P` (\eqn{p_{re} \times p_{re}} diagonal) from the validated `pfamily_list`: one \eqn{\tau^2_k} plug-in per component via `.two_block_tau2_plug_in_vector()`, assembled into `diag(tau2)` and inverted. `P` is no longer a caller-supplied argument on any `rLMMNormal_reg`/`rLMMindepNormalGamma_reg` export -- it is always this derived value, so it cannot be inconsistent with `pfamily_list`. |
| `.rLMM_validate_sim_method()` | `rLMM_reg.R` | Shared validator: `sim_method` must be exactly `"DEFAULT"` or `"TWO_BLOCK_GIBBS"`. Used by `rLMMNormal_reg_known_vcov()`, `rlmerb()`, and `rglmerb()`. |

### 4.4a The exact-iid engine (`sim_method = "DEFAULT"`)

| Function | File | Role |
|---|---|---|
| **`rLMMNormal_reg_known_vcov_iid()`** | `rLMM_reg.R` | Same matrix-input/dispersion/`pfamily_list` validation as `_two_bg()`, requires all-`dNormal`, then delegates to `.rLMMNormal_reg_run_iid()`. Unlike `_two_bg()`, it does **not** derive `P` itself -- `rLMMNormal_joint_iid()` below does that. |
| `.rLMMNormal_reg_run_iid()` | `rLMM_reg.R` | Builds a `P`-free `prior_list_block1 = list(dispersion, ddef = FALSE)`, calls `rLMMNormal_joint_iid()`, reshapes its output through the *same* `.rLMM_format_v2_out()` used by the Gibbs path (so downstream `fixef.*`/`coefficients`/`ranef.mode` staging is identical either way), then sets `m_convergence <- 1L`, `convergence_info$method <- "exact_iid"`, `draw_engine <- "rLMMNormal_joint_iid"`, and `sim_method_used <- "DEFAULT"`. |
| **`rLMMNormal_joint_iid()`** | `rLMMNormal_joint_iid.R` | Matrix-level export with the same `y, x, group, x_hyper, pfamily_list, prior_list_block1` signature as `two_block_rNormal_reg()`. Re-validates matrix inputs (`colnames(x)`, factor `group`) itself, then derives and injects `P` into `prior_list_block1` via `.rLMM_P_from_pfamily_list()` (`stop()`s if the caller already supplied `P`/`Sigma`) before `.two_block_validate_block1_prior()`. Builds the same `design`/`measurement_prior_list` shapes `.two_block_icm_at_start()` uses, then (1) `.lmerb_posterior_normal_system(design, measurement_prior_list)` builds `M` (posterior precision of the stacked Block~2 `gamma_full`) and, per group, `post_P_j` (conditional precision of `b_j \mid \gamma`), both independent of `gamma` (\S5); (2) `.lmerb_posterior_system_cholesky(system)` Cholesky-factors `M` and every `post_P_j` **once** (symmetry precondition below); (3) for each of the `n` draws, one `backsolve()` against `chol(M)` gives `gamma_full`, and one `backsolve()` per group against `chol(post_P_j)` gives `b_j` -- no iteration, no burn-in, no autocorrelation between draws. Returns `fixef_mean` (the same posterior mean `lmerb_posterior_mean()` would return) plus draw arrays in the shape `.two_block_format_cpp_out()` produces, so `.rLMM_format_v2_out()` (used by *both* engines) needs no engine-specific branching. |
| `.lmerb_posterior_normal_system()` | `lmebayes_posterior_icm.R` | Shared with `lmerb_posterior_mean()` (\S4.5) -- extracted from it without changing that function's signature/return. Same per-group `sigma2_j` indexing (\S5), so `"fixed_vector"` works identically on both engines. |
| `.lmerb_posterior_b_given_gamma()` | `lmebayes_posterior_icm.R` | Back-substitution step (mean of `b_j \mid \gamma`), also shared with `lmerb_posterior_mean()`; not needed for iid *sampling* itself (which draws `b_j` directly, \S4.4a step 3) but reused by `lmerb_posterior_mean()` for the posterior *mean*. |
| **`.lmerb_posterior_system_cholesky()`** | `lmebayes_posterior_icm.R` | Cholesky-factors `M` and every `post_P_j` from `.lmerb_posterior_normal_system()`'s output. **Validates `M` is numerically symmetric first** (within a relative tolerance) and defensively symmetrizes before factoring; `stop()`s with a message pointing at `sim_method = "TWO_BLOCK_GIBBS"` if the asymmetry is too large to safely ignore. `M` is only guaranteed exactly symmetric when `Sigma_ranef` (from `prior_list_block1`) is diagonal *and* its \eqn{k}-th diagonal entry equals component \eqn{k}'s `dNormal()` dispersion (\eqn{\tau^2_k}) in `pfamily_list` -- exactly how `lmerb()`/`glmerb()` always construct it, but not a precondition enforced at this matrix level, hence the runtime check. |

### 4.5 Shared run body: exact ICM + exact rate + sampling

| Function | File | Role |
|---|---|---|
| `.rLMMNormal_reg_run()` | `rLMM_reg.R` | Shared body for the `known_vcov` route (`estimated_vcov` uses the pilot-augmented `.rLMMNormal_reg_run_with_pilot()` instead). Builds `prior_list_block1` (with `P`, needed by `.rLMM_icm_at_start()`/`.rLM_calibrate_m_convergence()`), validates it, computes the ICM start, calibrates `m_convergence`, then calls `two_block_rNormal_reg()` with a **`P`-free** copy of `prior_list_block1` (that export derives its own `P`), stages output. |
| `.two_block_validate_block1_prior()` | `two_block_rNormal_reg.R` | Confirms `prior_list_block1` has `P`/`Sigma` and (for `gaussian()`) `dispersion`. Callers of the 5 matrix-level exports (`rLMMNormal_reg*`, `rGLMM_reg*`, `two_block_rNormal_reg()`, `rLMMNormal_joint_iid()`) must not pass `P`/`Sigma` themselves -- each export injects its own derived `P` right before calling this validator; internal-only callers (`.rLMM_icm_at_start()`, `.rLM_calibrate_m_convergence()`) still pass it directly. |
| `.rLMM_icm_at_start()` | `rLMM_reg.R` | Thin wrapper: builds the `design` list (`y, Z, groups, X_hyper, re_coef_names, group_name`) and calls `.two_block_icm_at_start()`; prints the "ICM posterior mean" verbose line. |
| `.two_block_icm_at_start()` | `two_block_measurement_prior.R` | Builds the `measurement_prior_list` (\S4.6) then, for `family = gaussian()`, calls **`lmerb_posterior_mean()`** (exact) instead of `glmerb_posterior_mode()` (ICM iteration, used for non-Gaussian families). |
| `.two_block_measurement_prior_list()` | `two_block_measurement_prior.R` | Converts `prior_list_block1$P`/`Sigma` + `dispersion` and each Block~2 `pfamily_list[[k]]$prior_list` into `list(dispersion_ranef, Sigma_ranef, prior_list)` -- `dispersion_ranef` here is exactly `prior_list_block1$dispersion` (scalar or named `J`-vector), passed straight through with `as.numeric()`. |
| **`lmerb_posterior_mean()`** | `lmebayes_posterior_icm.R` | **The exact multivariate-normal solver.** For each group `j`, forms `ZtZ_scaled[[j]] = crossprod(Z_j) / sigma2_j`, `Zty_scaled[[j]] = crossprod(Z_j, y_j) / sigma2_j`, where `sigma2_j <- if (length(sigma2) > 1) sigma2[[j]] else sigma2` (per-group value for `"fixed_vector"`, broadcast scalar for `"fixed"`). Eliminates all `b_j` algebraically (Schur complement; each `b_j` is affine in the shared Block~2 hyperparameter vector `gamma`, and never couples to another group's `b_{j'}`), leaving one linear system in `gamma` alone -- cost `O(J)` to build, no `J x J` or `J*p_re` matrix ever formed. Returns `fixef` (posterior mean of `gamma`), `b_mean` (posterior mean of every `b_j`), and `converged = TRUE, iterations = 1L, delta = 0` (no iteration occurs). |
| `build_mu_all()` | `build_mu_all.R` | Builds the group-by-group `mu_j(gamma)` design map (`X_hyper[[k]] %*% gamma_k` stacked per group) used inside `lmerb_posterior_mean()`'s Schur-complement algebra. Has a pure-R and a C++ (`use_cpp = TRUE` default) implementation. |
| `glmerb_posterior_mode()` | `lmebayes_posterior_icm.R` | *Not* on this route for Gaussian models (used for non-Gaussian GLMM families' ICM instead); mentioned because it shares `.two_block_measurement_prior_list()`/the same "shared Block~2 update" math and was fixed for the same per-group `sigma2_per_group[[j]]` indexing when Gaussian. |

### 4.6 Exact Theorem~3 rate calibration

| Function | File | Role |
|---|---|---|
| `.rLMM_calibrate_m_convergence()` | `rLMM_reg.R` | Calls `two_block_rate_from_pfamily_list()`, then `two_block_l_for_tv(rate, tv_tol, method = "theorem3")` + `.two_block_cap_inner_sweeps()` to get `m_convergence`; tags `convergence_info$method` via `.rLMM_rate_calibration_meta()` -- `"exact"` when `any_non_normal = FALSE` and no random measurement dispersion is involved (always true on this route: fixed dispersion + all-dNormal). |
| `two_block_rate_from_pfamily_list()` | `two_block_ergodicity.R` | Thin adapter from the `pfamily_list`/`prior_list_block1` shapes to `two_block_rate()`'s raw-matrix interface. |
| `two_block_rate()` | `two_block_ergodicity.R` | Orchestrates the eigenvalue computation: `.two_block_rate_inputs()` -> `.two_block_S_P11()` -> `.two_block_gen_eigen()`; returns `lambda_star` (spectral radius / convergence rate) and `eigenvalues`. |
| `.two_block_rate_inputs()` | `two_block_ergodicity.R` | Builds per-observation working weights from `prior_list_block1$dispersion`: if `length(disp) == J`, expands the per-group vector to one weight per observation via `block_info$rows`/`row_idx` (the bug fixed for `"fixed_vector"`); a length-1 `disp` still broadcasts to every observation as before. |
| `.two_block_S_P11()` | `two_block_ergodicity.R` | Forms the Schur-complement-style matrices needed for the generalized eigenproblem. |
| `.two_block_gen_eigen()` | `two_block_ergodicity.R` | Solves the generalized eigenvalue problem; the largest eigenvalue is the two-block Gibbs sampler's asymptotic convergence rate `lambda_star`. |
| `two_block_l_for_tv()` | `two_block_ergodicity.R` | Converts `lambda_star` + `tv_tol` into a minimum sweep count via Theorem~3's total-variation bound. |
| `.two_block_cap_inner_sweeps()` | `two_block_ergodicity.R` | Applies a sane upper cap to the computed `m_min` before it becomes `m_convergence`. |

### 4.7 Sampling orchestration and output staging

| Function | File | Role |
|---|---|---|
| `two_block_rNormal_reg()` | `two_block_rNormal_reg.R` | Runs `n` replicate draws, each a full pass of `m_convergence` inner two-block Gibbs sweeps starting from the exact ICM mean; delegates the actual sweep loop to compiled code via `.two_block_rNormal_reg_cpp()`. |
| `.two_block_rNormal_reg_cpp()` | `rcpp_wrappers.R` | R-to-C++ bridge: marshals `y, x/Z, block, x_hyper, prior_list_block1, pfamily_list, fixef_start, m_convergence, ...` into the arguments expected by the compiled kernel (\S6) and marshals the result back into R lists/matrices. |
| `.rLMM_format_v2_out()` | `rLMM_reg.R` | Reshapes the raw sampler output into the staged `fixef.*` / `coefficients` / `ranef.mode` namespaces expected by `rlmerb()`/`rglmerb()`. |
| `.two_block_as_staged_names()` | `two_block_glmm_pilot_helpers.R` | Name-staging helper shared by both the fixed and pilot output paths. |
| `.lmebayes_add_fixef_summaries()` | `mixed_rmerb_helpers.R` | Adds posterior summary statistics (mean/sd/quantiles) for each Block~2 fixed-effect vector to the returned object. |

---

## 5. Why this route is exact (math sketch)

With fixed \eqn{\sigma^2} (scalar or per-group) and all-`dNormal` Block~2
priors, the model is:

\deqn{
  y_j \mid b_j \sim N(Z_j b_j,\, \sigma^2_j I), \qquad
  b_j \mid \gamma \sim N(\mu_j(\gamma),\, \Sigma_b), \qquad
  \gamma_k \sim N(\mu_{\gamma_k},\, \Sigma_{\gamma_k})
}

for group \eqn{j} and Block~2 component \eqn{k}, with \eqn{\mu_j(\gamma)}
linear in the stacked \eqn{\gamma} (via `build_mu_all()`'s
`X_hyper[[k]] %*% gamma_k` map). Every density here is Gaussian, so the joint
`(gamma, b_1, ..., b_J)` posterior is exactly multivariate normal:

- **Posterior mean = posterior mode**, computable in one linear solve
  (`lmerb_posterior_mean()`, \S4.5) -- no ICM *iteration* is needed, unlike
  the non-Gaussian (`glmerb_posterior_mode()`) or estimated-\eqn{\tau^2}
  routes, which must alternate Block~1/Block~2 updates to a tolerance.
- **The two-block Gibbs sampler's convergence rate is exact**, not a
  conservative bound -- Theorem~3's rate depends only on the *linear*
  Gibbs-sweep operator built from `P_b`, `P_{\gamma_k}`, `sigma2_j`, and the
  design (`.two_block_gen_eigen()`, \S4.6); there is no nonlinearity (GLM
  link, non-Gaussian RE prior) to bound conservatively.

`"fixed_vector"` differs from `"fixed"` only in letting `sigma2_j` vary by
group instead of forcing `sigma2_j = sigma2` for all `j`; every formula above
already carries the group subscript, so no new machinery was needed --
only the `sigma2_j` group-indexing itself, which is what the two ICM/rate
fixes (\S2) corrected.

---

## 6. C++ kernel chain (`src/`)

`two_block_rNormal_reg()`'s inner sweep loop (Block~1 draws, `is_ing = false`
Block~2 draws) is compiled; the fixed/known-vcov route reaches:

```
.two_block_rNormal_reg_cpp()                 [rcpp_wrappers.R -- R entry point]
  -> block1_prior_list()                     [twoBlockGibbs.cpp -- Block~1 prior payload]
  -> MuAllBuilder::build()                   [simfuncs.h, impl in twoBlockGibbs.cpp -- mu_j(gamma) map]
  -> block_rNormalReg_cpp_export()           [block_utils.cpp -- per-sweep Block~1 draw]
       -> rNormalRegBlocks()                 [rNormalRegBlocks.cpp -- loops groups]
            -> rNormalReg()                  [rNormalReg.cpp -- one group's Gaussian draw, no lm.fit]
  -> (Block~2 loop, is_ing = false branch)
       -> rNormalReg()                       [rNormalReg.cpp -- dNormal conjugate draw per component]
```

No accept/reject envelope, no `glmbayesCore::rglmb()` call, and no
`lm.fit()` appear anywhere in this chain (see
[`DESIGN_RNORMALREG_BLOCKS.md`](DESIGN_RNORMALREG_BLOCKS.md) for the
Gaussian-block-Gibbs migration that removed the latter). Every draw is a
direct conjugate-normal sample.

---

## 7. Related files

| Topic | Path |
|---|---|
| Estimated-\eqn{\tau^2} sibling route (pilot stage, `_estimated_vcov`) | `R/rLMM_reg.R` (`rLMMNormal_reg_estimated_vcov`, `.rLMMNormal_reg_run_with_pilot`) |
| Sampled-dispersion routes (`"gamma"`/`"gamma_list"`) | `R/rLMM_reg.R` (`rLMMindepNormalGamma_reg_*`), `inst/BLOCK_ING_RINDEPNORMALGAMMA_REG.md` |
| Convergence theory background | `inst/BLOCK_GIBBS_ERGODICITY.md` |
| Gaussian block-Gibbs C++ migration (Block~1 kernel) | `inst/DESIGN_RNORMALREG_BLOCKS.md` |
| Full export inventory / reachability | `inst/R_EXPORT_REACHABILITY.md`, `inst/R_FUNCTION_INVENTORY.md` |
| `"fixed_vector"` feature history | `NEWS.md` (this package and `lmebayes`) |
| `sim_method` / exact-iid engine feature history | `NEWS.md` (this package and `lmebayes`) |
| Exact-iid sampler implementation | `R/rLMMNormal_joint_iid.R` (`rLMMNormal_joint_iid`), `R/lmebayes_posterior_icm.R` (`.lmerb_posterior_normal_system`, `.lmerb_posterior_b_given_gamma`, `.lmerb_posterior_system_cholesky`) |

---

## 8. Verification method

- Read `lmebayesCore/R/rLMM_reg.R` in full for `rLMMNormal_reg`,
  `rLMMNormal_reg_known_vcov`, `.rLMMNormal_reg_run`, `.rLMM_icm_at_start`,
  `.rLMM_calibrate_m_convergence`, and every `.rLMM_validate_*` helper it
  calls.
- Read `R/two_block_measurement_prior.R`
  (`.two_block_measurement_prior_list`, `.two_block_icm_at_start`) and
  `R/lmebayes_posterior_icm.R` (`lmerb_posterior_mean`,
  `glmerb_posterior_mode`) in full, including the per-group `sigma2_j`
  indexing fixed for `"fixed_vector"`.
- Read `R/two_block_ergodicity.R`
  (`.two_block_rate_inputs`, `.two_block_S_P11`, `.two_block_gen_eigen`,
  `two_block_rate`, `two_block_rate_from_pfamily_list`,
  `two_block_l_for_tv`, `.two_block_cap_inner_sweeps`), including the
  per-group weight expansion fixed for `"fixed_vector"`.
- Read `R/two_block_rNormal_reg.R`
  (`two_block_rNormal_reg`, `.two_block_validate_block1_prior`,
  `.two_block_validate_pfamily_list`, `.two_block_summarize_pfamily_list`)
  and `R/rcpp_wrappers.R` (`.two_block_rNormal_reg_cpp`).
- Read `R/mixed_rmerb_helpers.R` in full for the dispersion resolver family
  (`.lmebayes_resolve_dispersion_ranef*`), the route-glue helpers
  (`.lmebayes_matrix_args_lmm`, `.lmebayes_run_lmm_engine`,
  `.lmebayes_block1_prior_list`, `.lmebayes_attach_sigma2`), and
  `R/lmebayes_reg_route_table.R` (`.lmebayes_reg_route_key`,
  `REG_ROUTE_TABLE`, `.lmebayes_reg_route_fn`).
- Read `R/rlmerb.R` / `R/rglmerb.R` end to end for the matrix-level entry
  points, including the `sim_method` formal and `.rLMM_validate_sim_method()`.
- Read `R/rLMMNormal_joint_iid.R` (`rLMMNormal_joint_iid`) end to end, and the
  `rLMMNormal_reg_known_vcov()`/`_iid()`/`_two_bg()`/`.rLMMNormal_reg_run_iid()`
  split in `R/rLMM_reg.R`, cross-checking `.lmerb_posterior_normal_system()`/
  `.lmerb_posterior_system_cholesky()` in `R/lmebayes_posterior_icm.R` for the
  `M` symmetry precondition.
- Grepped `src/*.cpp` for `block1_prior_list`, `MuAllBuilder`,
  `block_rNormalReg_cpp_export`, `rNormalRegBlocks`, `rNormalReg` to confirm
  the compiled kernel chain names and files.

---

## 9. Changelog

| Date | Note |
|---|---|
| 2026-07-20 | Initial doc: traced the full `known_vcov` (fixed dispersion + all-`dNormal` Block~2) call graph from `lmerb()`/`glmerb()` through `rLMMNormal_reg_known_vcov()`'s exact ICM (`lmerb_posterior_mean()`) and exact Theorem~3 rate calibration, down to the compiled `rNormalReg`/`rNormalRegBlocks` kernels. Written alongside the new `"fixed_vector"` (per-group fixed dispersion) `dispersion_ranef` mode, which reuses this exact route unchanged apart from `sigma2_j` becoming per-group. |
| 2026-07-20 | Added \S4.4a: the new `sim_method = "DEFAULT"` exact-iid engine (`rLMMNormal_reg_known_vcov_iid()` / `rLMMNormal_joint_iid()`), which draws directly from the same closed-form Gaussian that `lmerb_posterior_mean()` already summarizes exactly (no Gibbs sweeps). The pre-existing Gibbs body is now `rLMMNormal_reg_known_vcov_two_bg()`; `rLMMNormal_reg_known_vcov()` is a thin `sim_method` dispatcher between the two. Updated \S1-\S3, \S4.1, \S4.3, \S4.4, \S7, \S8 accordingly. |
