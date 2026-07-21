# Front-door functions: `lmerb()` / `glmerb()` shared entry path

**Scope.** This document enumerates the functions -- R user-facing and R
internal (`@noRd`) -- that form the **shared entry path** used by
`lmerb()` and `glmerb()` (package **lmebayes**) regardless of which
`dispersion_ranef` mode or variance-component mode the model ends up in: it
covers everything from the formula call down to the point where the flat,
matrix-level arguments are handed to a specific route export (e.g.
`rLMMNormal_reg_known_vcov()`).

This is deliberately **not** specific to any one route. In particular:

- `model_setup()`, `Prior_Setup_lmebayes()`, `pfamily_list()`,
  `.lmebayes_priors_from_pfamily_list()`, and the `dispersion_ranef`
  resolution helpers run identically for `disp_mode` in
  `{"fixed", "fixed_vector", "gamma", "gamma_list"}`.
- `.lmebayes_reg_route_key()` / `.lmebayes_reg_route_fn()` (via
  `REG_ROUTE_TABLE`, `R/lmebayes_reg_route_table.R`) resolve to **any** of
  the six registered routes (`lmm_fixed_known`, `lmm_fixed_estimated`,
  `lmm_gamma_known`, `lmm_gamma_estimated`, `glmm_known`, `glmm_estimated`),
  not only the fixed-dispersion + known-vcov one.
- `.lmebayes_run_lmm_engine()` (and its arg-builder,
  `.lmebayes_matrix_args_lmm()`) is shared by all four **LMM** routes (fixed
  or gamma dispersion, known or estimated vcov); the parallel
  `.lmebayes_run_glmm_engine()` / `.lmebayes_matrix_args_glmm()` pair (not
  detailed here) handles the two non-Gaussian **GLMM** routes.

The route-specific engine reached when `disp_mode` is fixed (scalar or
per-group vector) **and** every `pfamily_list` component is `dNormal()` --
i.e. `rLMMNormal_reg_known_vcov()` and everything downstream of it (both the
exact-iid and two-block-Gibbs sub-engines) -- is documented separately in
`inst/README_LMM_FIXED_DISPERSION_KNOWN_VCOV_FUNCTIONS.md`. The other five
routes are not documented in either file as of this writing.

Every function name and file path below was confirmed by reading the current
source (not inferred).

---

## 1. Call chain

```
lmerb() / glmerb()                              [lmebayes]
  -> model_setup()                               [lmebayesCore]  (design: y, Z, groups, X_hyper, ranks; optional lmer/glmer ref fit)
       -> .lmebayes_mer_optional_args()          [lmebayesCore]  (optional args forwarded to lmer()/glmer() for the reference fit)
       -> extract_re_hyper_matrices()            [lmebayesCore]  (Z / X_hyper from the reference fit)
       -> extract_mer_variance_components()       [lmebayesCore]  (variance components from the reference fit)
  -> Prior_Setup_lmebayes()                      [lmebayesCore]  (calibrate Block-1 Sigma_ranef / dispersion, Block-2 hyperpriors)
  -> pfamily_list(<lmebayes_prior_setup>)        [lmebayesCore]  (builds glmbayesCore::dNormal(...) objects)
  -> lmebayesCore::rlmerb() / rglmerb()
       -> .lmebayes_priors_from_pfamily_list()
            -> .lmebayes_resolve_dispersion_ranef()
                 -> .lmebayes_resolve_dispersion_ranef_fixed_vector()   (named-vector dispersion_ranef only)
       -> .lmebayes_run_lmm_engine()              (Gaussian families only; non-Gaussian uses .lmebayes_run_glmm_engine(), not covered here)
            -> .lmebayes_reg_route_key()
            -> .lmebayes_reg_route_fn()           -> resolves REG_ROUTE_TABLE[[route_key]]$export to a callable
            -> .lmebayes_matrix_args_lmm()
            -> do.call(<route export>, args)      <-- e.g. rLMMNormal_reg_known_vcov() for the fixed+known route
            -> .lmebayes_attach_sigma2()          (post-processes the route's raw output)
       -> [rlmerb()/rglmerb() only, if verbose]
            .lmebayes_block2_icm_labels() -> .lmebayes_print_icm_fixef_table()
       -> .lmebayes_add_fixef_summaries()         (adds fixef.means etc. to the final output)
```

`glmerb(family = gaussian())` does not take a separate GLMM code path: for
the Gaussian family it is routed through the identical
`.lmebayes_run_lmm_engine()` call used by `rlmerb()`.

---

## 2. R user-facing (exported) functions

| Package | Function | File | Role |
|---|---|---|---|
| lmebayes | `lmerb()` | `R/lmerb.R` | Formula entry point; builds design/prior, calls `rlmerb()`. |
| lmebayes | `glmerb()` | `R/glmerb.R` | Formula entry point; for `family = gaussian()` routes to the same LMM engine as `lmerb()`; otherwise the (not covered here) GLMM engine. |
| lmebayesCore | `model_setup()` | `R/model_setup.R` | Builds `design` (y, Z, groups, X_hyper, ranks) from a formula; optional `lmer`/`glmer` reference fit. |
| lmebayesCore | `Prior_Setup_lmebayes()` | `R/Prior_Setup_lmebayes.R` | Calibrates Block-1 `Sigma_ranef`/dispersion and Block-2 hyperpriors from the reference fit. |
| lmebayesCore | `pfamily_list()` (S3, `lmebayes_prior_setup` method) | `R/pfamily_list_lmebayes_prior_setup.R` | Converts calibrated priors into a list of `glmbayesCore::dNormal()` objects. |
| lmebayesCore | `rlmerb()` | `R/rlmerb.R` | Matrix-level LMM entry; validates inputs, calls `.lmebayes_run_lmm_engine()`, applies verbose ICM printing and fixef summaries. |
| lmebayesCore | `rglmerb()` | `R/rglmerb.R` | Matrix-level GLMM entry; for Gaussian family calls the identical `.lmebayes_run_lmm_engine()`. |
| lmebayesCore | `extract_mer_variance_components()` | `R/lme4_design_utilities.R` | Extracts variance components from an `lmer`/`glmer` reference fit (used by `model_setup()`). |
| lmebayesCore | `extract_re_hyper_matrices()` | `R/lme4_design_utilities.R` | Extracts `Z`/`X_hyper` design pieces from an `lmer`/`glmer` reference fit (used by `model_setup()`). |
| glmbayesCore | `dNormal()` | (glmbayesCore) | Constructs the Block-2 prior objects built by `pfamily_list()`. |

## 3. R internal (`@noRd`) helpers

| Function | File | Role | Calls |
|---|---|---|---|
| `.lmebayes_mer_optional_args()` | `R/model_setup.R` | Builds optional-argument list forwarded to `lmer()`/`glmer()` for the reference fit inside `model_setup()`. | -- |
| `.lmebayes_resolve_dispersion_ranef()` | `R/mixed_rmerb_helpers.R` | Dispatches `dispersion_ranef` to one of `none`/`fixed`/`fixed_vector`/`gamma`/`gamma_list` modes. | `.lmebayes_resolve_dispersion_ranef_fixed_vector()` (named-vector case) |
| `.lmebayes_resolve_dispersion_ranef_fixed_vector()` | `R/mixed_rmerb_helpers.R` | Validates/aligns a named per-group fixed-dispersion vector against `design$groups` levels. | -- |
| `.lmebayes_priors_from_pfamily_list()` | `R/mixed_rmerb_helpers.R` | Normalizes `pfamily_list` + `dispersion_ranef` into the `prior` object (`Sigma_ranef`, `prior_list`, `ptypes`, `any_non_normal`, ...). | `.lmebayes_resolve_dispersion_ranef()` |
| `.lmebayes_matrix_args_lmm()` | `R/mixed_rmerb_helpers.R` | Assembles the flat argument list passed to the routed export (incl. setting `sim_method` only for the `lmm_fixed_known` route). | -- |
| `.lmebayes_run_lmm_engine()` | `R/mixed_rmerb_helpers.R` | Resolves the route key, builds args, `do.call()`s the route export, attaches `sigma2`/`sigma2.mean` and a default `sim_method_used`. | `.lmebayes_reg_route_key()`, `.lmebayes_reg_route_fn()`, `.lmebayes_matrix_args_lmm()`, `.lmebayes_attach_sigma2()` |
| `.lmebayes_attach_sigma2()` | `R/mixed_rmerb_helpers.R` | Post-processes sampler output to attach `sigma2`/`sigma2.mean` per dispersion mode (`fixed`/`fixed_vector` are pass-through constants). | -- |
| `.lmebayes_add_fixef_summaries()` | `R/mixed_rmerb_helpers.R` | Adds posterior-mean summaries (`fixef.means`, etc.) to the final `rlmerb()`/`rglmerb()` output. | -- |
| `.lmebayes_block2_icm_labels()` / `.lmebayes_print_icm_fixef_table()` | `R/mixed_rmerb_helpers.R` | Verbose-mode ICM reporting labels/table, printed from `rlmerb()`/`rglmerb()`. | -- |
| `.lmebayes_reg_route_key()` | `R/lmebayes_reg_route_table.R` | Maps `(family, disp_mode, any_non_normal)` -> one of the six route-key strings. | -- |
| `.lmebayes_reg_route_fn()` | `R/lmebayes_reg_route_table.R` | Resolves a route key to `{export_fn, needs_pilot, draw_engine_label}` via the `REG_ROUTE_TABLE` list. | -- |

`.lmebayes_run_glmm_engine()` / `.lmebayes_matrix_args_glmm()` /
`.lmebayes_block1_prior_list()` (also in `R/mixed_rmerb_helpers.R`) are the
GLMM (non-Gaussian) counterparts of `.lmebayes_run_lmm_engine()` /
`.lmebayes_matrix_args_lmm()`; they are not detailed here since no route they
resolve to is documented in this pair of files.

---

## 4. Verification method / caveats

Every function/file pairing above was confirmed against the source in this
workspace (`c:\Rpackages\lmebayesCore`, `c:\Rpackages\lmebayes`) by direct
reading and `grep`, not inferred from documentation or memory. In particular,
the `rlmerb()`/`rglmerb()` call sites of `.lmebayes_block2_icm_labels()`,
`.lmebayes_print_icm_fixef_table()`, and `.lmebayes_add_fixef_summaries()`
were confirmed by reading `R/rlmerb.R` and `R/rglmerb.R` directly.
