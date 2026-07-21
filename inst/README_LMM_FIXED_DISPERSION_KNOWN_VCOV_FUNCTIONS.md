# Function inventory: LMM models with fixed dispersion and known variance components

**Scope.** This document enumerates the functions -- R user-facing, R internal
(`@noRd`), and C++ -- that are **specific to** the `lmerb()` /
`glmerb(family = gaussian())` route for models with:

- `dispersion_ranef` **fixed** (a positive scalar, or a named per-group
  numeric vector -- `disp_mode == "fixed"` / `"fixed_vector"`), and
- every component of `pfamily_list` a `dNormal()` (no
  `dIndependent_Normal_Gamma` components -- **known**, not estimated,
  variance components).

This is the route internally keyed `lmm_fixed_known` (see `REG_ROUTE_TABLE` in
`R/lmebayes_reg_route_table.R`). The statistical description of this route
lives in `inst/README_KNOWN_VCOV_GAUSSIAN.md`; this document is the matching
**engineering inventory**, for use if/when this model's functions are ported
into **glmbayesCore**. Every function name and file path below was confirmed
by reading the current source (not inferred).

**Out of scope for this document** (tracked separately, not part of "this
model"):

- The shared entry path used by `lmerb()`/`glmerb()` for **every** route
  (formula/design setup, prior calibration, `dispersion_ranef` resolution,
  route-table lookup, and matrix-argument assembly) -- documented separately
  in `inst/README_LMERB_GLMERB_FRONT_DOOR.md`. This document picks up at
  `rLMMNormal_reg_known_vcov()`, the point where that shared front door hands
  off to the route specific to this model.
- The generic per-group/per-block conjugate Normal-regression sampling
  primitives this route's Gibbs engine calls into (`rNormalReg()`,
  `rNormalRegBlocks()`, `block_rNormalReg_cpp_export()`, `normalize_block()`).
  These are shared infrastructure also used by (or duplicative of) other
  parts of the `glmbayes`/`glmbayesCore` ecosystem (e.g. the machinery behind
  `lmbBlock()`/`glmbBlock()`); whether/how to consolidate them is a separate,
  general "block-sampling primitives" question, independent of this model.
  They are listed once in Section 3 as dependencies, with no further
  porting discussion here.
- Any other function or route not reached by this specific model (ING
  routes, GLMM non-Gaussian routes, envelope/OpenCL machinery, etc.).

The route has two sub-engines, selected by the `sim_method` argument on
`lmerb()`/`glmerb()`/`rlmerb()`/`rglmerb()`:

- `sim_method = "DEFAULT"` -- exact **iid** draws from the closed-form joint
  Gaussian posterior (no Markov chain).
- `sim_method = "TWO_BLOCK_GIBBS"` -- the two-block Gibbs sampler.

Both are reachable from the same dispatcher, `rLMMNormal_reg_known_vcov()`.

---

## 1. Call chains

Both sub-engines below are reached the same way: `lmerb()`/`glmerb()` ->
... -> `.lmebayes_run_lmm_engine()` -> `rLMMNormal_reg_known_vcov()`. That
shared front door (`model_setup()`, `Prior_Setup_lmebayes()`,
`pfamily_list()`, `.lmebayes_priors_from_pfamily_list()`,
`.lmebayes_run_lmm_engine()`, route-table lookup, `.lmebayes_matrix_args_lmm()`,
etc.) is documented in `inst/README_LMERB_GLMERB_FRONT_DOOR.md`; the chains
below start at `rLMMNormal_reg_known_vcov()`, the function that document hands
off to for this model.

### 1.1 `sim_method = "DEFAULT"` -- exact iid engine

```
rLMMNormal_reg_known_vcov()
  -> rLMMNormal_reg_known_vcov_iid()
       -> .rLMMNormal_reg_run_iid()
            -> rLMMNormal_joint_iid()
                 -> .lmerb_posterior_normal_system()     (pure R: solve/crossprod -- builds M, per-group post_P_j, v)
                 -> .lmerb_posterior_system_cholesky()   (pure R: chol(), with symmetry guard)
                 -> [per draw] rnorm() + backsolve() for gamma_k; build_mu_all(); backsolve() for each b_j
            -> .rLMM_format_v2_out() -> .two_block_as_staged_names()
```

All of the joint-Gaussian sampling math (`solve`, `chol`, `backsolve`,
`rnorm`) is base R and specific to this model. The only C++ touched is
`build_mu_all()`'s default backend (`MuAllBuilder`, Section 3) -- a shared
dependency, not part of this model's own code.

### 1.2 `sim_method = "TWO_BLOCK_GIBBS"` -- two-block Gibbs engine

```
rLMMNormal_reg_known_vcov()
  -> rLMMNormal_reg_known_vcov_two_bg()
       -> .rLMMNormal_reg_run()
            -> .rLMM_icm_at_start() -> .two_block_icm_at_start()
                 -> .two_block_measurement_prior_list()
                 -> lmerb_posterior_mean()
                      -> .lmerb_posterior_normal_system()
                      -> .lmerb_posterior_b_given_gamma()   -> build_mu_all()  [shared dependency, Section 3]
            -> .rLMM_calibrate_m_convergence() -> .rLMM_rate_calibration_meta()
                 -> two_block_rate_from_pfamily_list() -> two_block_rate()
                      -> .two_block_rate_inputs()  (uses normalize_block() -- shared dependency, Section 3)
                      -> .two_block_S_P11() -> .two_block_gen_eigen()
                 -> two_block_l_for_tv() -> two_block_tv_bound() -> .two_block_tv_bound_one() / .two_block_erfn()
                 -> .two_block_cap_inner_sweeps()
            -> two_block_rNormal_reg()
                 -> glmbayesCore::glmbfamfunc()          (R-level only; f2/f3/f2_gauss/f3_gauss)
                 -> .two_block_rNormal_reg_cpp()
                      -> .Call(`_lmebayesCore_two_block_rNormal_reg_v2_cpp_export`)
            -> .rLMM_format_v2_out() -> .two_block_as_staged_names()
```

Inside `two_block_rNormal_reg_v2_cpp_export()` (`src/twoBlockGibbs.cpp`), each
outer draw runs `m_convergence` inner sweeps of Block 1 (random effects,
delegated to the shared `block_rNormalReg_cpp_export()` primitive) and Block 2
(hyper-parameters, delegated to the shared `rNormalReg()` primitive). The
model-specific parts of that driver are the two-block bookkeeping around those
calls (`MuAllBuilder`/`block1_prior_list()` rebuilding the Block-1 prior from
`mu_all` each sweep, and `block2_prior_prep_v2()` parsing the `dNormal`
pfamilies) -- the actual per-group conjugate-Normal draw is the shared
primitive described in Section 3.

---

## 2. Functions specific to this model

### 2.1 R user-facing (exported)

| Package | Function | File | Role |
|---|---|---|---|
| lmebayesCore | `rLMMNormal_reg()` | `R/rLMM_reg.R` | Higher dispatcher (known vs. estimated vcov); calls `rLMMNormal_reg_known_vcov()` when vcov is known. Not itself on the `lmerb()`/`glmerb()` path (which reaches `rLMMNormal_reg_known_vcov()` directly via the route table -- see `inst/README_LMERB_GLMERB_FRONT_DOOR.md`), but listed here since it shares this model's code. |
| lmebayesCore | `rLMMNormal_reg_known_vcov()` | `R/rLMM_reg.R` | **`sim_method` dispatcher** -- routes to `_iid` or `_two_bg`. This is the hand-off point from the shared front door (`inst/README_LMERB_GLMERB_FRONT_DOOR.md`) into this model. |
| lmebayesCore | `rLMMNormal_reg_known_vcov_iid()` | `R/rLMM_reg.R` | Exact-iid engine entry point. |
| lmebayesCore | `rLMMNormal_reg_known_vcov_two_bg()` | `R/rLMM_reg.R` | Two-block Gibbs engine entry point. |
| lmebayesCore | `rLMMNormal_joint_iid()` | `R/rLMMNormal_joint_iid.R` | Closed-form MVN construction + iid draws (algorithmic core of the iid engine). |
| lmebayesCore | `lmerb_posterior_mean()` | `R/lmebayes_posterior_icm.R` | Exact joint posterior mean of `(gamma, b)`; the closed-form answer for the iid engine and the ICM start for the Gibbs engine. |
| lmebayesCore | `two_block_rNormal_reg()` | `R/two_block_rNormal_reg.R` | R-level orchestration of the two-block Gibbs C++ driver (family callbacks, argument marshaling). |
| lmebayesCore | `two_block_rate()` | `R/two_block_ergodicity.R` | Theorem-3 convergence-rate eigenvalue computation from raw block1/block2 priors. |
| lmebayesCore | `two_block_rate_from_pfamily_list()` | `R/two_block_ergodicity.R` | Adapter: `pfamily_list`-shaped priors -> `two_block_rate()`. |
| lmebayesCore | `two_block_tv_bound()` | `R/two_block_ergodicity.R` | Total-variation distance bound from the calibrated rate at a given number of sweeps. |
| lmebayesCore | `two_block_l_for_tv()` | `R/two_block_ergodicity.R` | Minimum number of inner sweeps to hit a target `tv_tol`. |
| glmbayesCore | `glmbfamfunc()` | (glmbayesCore) | Supplies `f2`/`f3` family callbacks to the Gibbs C++ driver; marshaled through but not exercised by the Gaussian conjugate branch. |

`build_mu_all()` is used heavily by this model but is a shared dependency, not
specific to it -- see Section 3. `lmerb()`, `glmerb()`, `rlmerb()`,
`rglmerb()`, `model_setup()`, `Prior_Setup_lmebayes()`, `pfamily_list()`, and
`glmbayesCore::dNormal()` are the shared front door and are documented in
`inst/README_LMERB_GLMERB_FRONT_DOOR.md`, not here.

### 2.2 R internal (`@noRd`) helpers

| Function | File | Role | Calls (within this model) |
|---|---|---|---|
| `.rLMM_validate_sim_method()` | `R/rLMM_reg.R` | Validates `sim_method %in% c("DEFAULT", "TWO_BLOCK_GIBBS")`. | -- |
| `.rLMM_validate_matrix_inputs()` | `R/rLMM_reg.R` | Shared matrix/shape validation for `y`, `x`, `block`, `x_hyper`, etc. across `rLMMNormal_reg_*`. | -- |
| `.rLMMNormal_reg_run_iid()` | `R/rLMM_reg.R` | Body of the iid engine: calls `rLMMNormal_joint_iid()` and formats the result to match Gibbs-engine output shape. | `rLMMNormal_joint_iid()`, `.rLMM_format_v2_out()` |
| `.rLMMNormal_reg_run()` | `R/rLMM_reg.R` | Body of the Gibbs engine: ICM start, rate calibration, `two_block_rNormal_reg()`, output formatting. | `.rLMM_icm_at_start()`, `.rLMM_calibrate_m_convergence()`, `two_block_rNormal_reg()`, `.rLMM_format_v2_out()` |
| `.rLMM_icm_at_start()` | `R/rLMM_reg.R` | Thin wrapper calling `.two_block_icm_at_start()` for the Gaussian-known case. | `.two_block_icm_at_start()` |
| `.rLMM_calibrate_m_convergence()` / `.rLMM_rate_calibration_meta()` | `R/rLMM_reg.R` | Computes the exact Theorem-3 number of inner Gibbs sweeps (`m_convergence`) from `tv_tol`. | `two_block_rate_from_pfamily_list()`, `two_block_l_for_tv()`, `.two_block_cap_inner_sweeps()` |
| `.rLMM_format_v2_out()` | `R/rLMM_reg.R` | Normalizes either engine's raw output into the staged (`fixef.*`, `ranef.*`, ...) public shape. | `.two_block_as_staged_names()` |
| `.lmerb_posterior_normal_system()` | `R/lmebayes_posterior_icm.R` | Builds the joint linear system (`M`, per-group `post_P_j`, `v`) for the exact Gaussian posterior of `(gamma, b)`. | -- |
| `.lmerb_posterior_b_given_gamma()` | `R/lmebayes_posterior_icm.R` | Computes `E[b \| gamma]` given a `gamma` draw/estimate. | `build_mu_all()` (shared dependency) |
| `.lmerb_posterior_system_cholesky()` | `R/lmebayes_posterior_icm.R` | Cholesky factors of `M` and each `post_P_j`, with a symmetry validation/repair guard. | -- |
| `.two_block_validate_pfamily_list()` / `.two_block_summarize_pfamily_list()` | `R/two_block_rNormal_reg.R` | Validates `pfamily_list` structure; on this route, asserts the "all dNormal" gate. | -- |
| `.two_block_validate_block1_prior()` / `.two_block_normalize_family()` | `R/two_block_rNormal_reg.R` | Validates the Block-1 `{P, dispersion, ddef}` prior; normalizes the `family` argument. | -- |
| `.two_block_format_cpp_out()` | `R/two_block_rNormal_reg.R` | Converts the raw C++ list (`fixef_draws`, `b_arr`, ...) into the R-facing list shape. | -- |
| `.two_block_rNormal_reg_cpp()` | `R/rcpp_wrappers.R` | Thin positional `.Call()` bridge to `_lmebayesCore_two_block_rNormal_reg_v2_cpp_export`. | C++ |
| `.two_block_measurement_prior_list()` | `R/two_block_measurement_prior.R` | Builds the Block-1 measurement prior list (`Sigma_ranef`, dispersion) consumed by the ICM start and the Gibbs engine. | -- |
| `.two_block_icm_at_start()` | `R/two_block_measurement_prior.R` | ICM-at-start orchestration for the Gaussian-known case; calls `lmerb_posterior_mean()` directly (no iteration needed -- exact mean). | `.two_block_measurement_prior_list()`, `lmerb_posterior_mean()` |
| `.two_block_rate_inputs()` | `R/two_block_ergodicity.R` | Normalizes raw block1/block2 prior inputs into the matrices consumed by `.two_block_S_P11()`/`.two_block_gen_eigen()`. | `normalize_block()` (shared dependency) |
| `.two_block_S_P11()` | `R/two_block_ergodicity.R` | Assembles the `S` (likelihood) and `P11` (Block-1 prior precision) matrices for the rate eigenproblem. | -- |
| `.two_block_gen_eigen()` | `R/two_block_ergodicity.R` | Generalized eigenvalue computation underlying the Theorem-3 rate. | -- |
| `.two_block_tv_bound_one()` / `.two_block_erfn()` | `R/two_block_ergodicity.R` | Closed-form single-sweep TV-bound helpers. | -- |
| `.two_block_cap_inner_sweeps()` | `R/two_block_ergodicity.R` | Coerces/caps the calibrated sweep count to a valid integer range. | -- |
| `.two_block_uncertified_l_fallback()` | `R/two_block_ergodicity.R` | Fallback sweep-count heuristic when the certified bound is not applicable. | -- |
| `.two_block_as_staged_names()` | `R/two_block_glmm_pilot_helpers.R` | Final output-name staging shared by both engines' formatting step. | -- |

`.lmebayes_resolve_dispersion_ranef()`, `.lmebayes_resolve_dispersion_ranef_fixed_vector()`,
`.lmebayes_priors_from_pfamily_list()`, `.lmebayes_matrix_args_lmm()`,
`.lmebayes_run_lmm_engine()`, `.lmebayes_attach_sigma2()`,
`.lmebayes_add_fixef_summaries()`, `.lmebayes_block2_icm_labels()`,
`.lmebayes_print_icm_fixef_table()`, `.lmebayes_reg_route_key()`,
`.lmebayes_reg_route_fn()`, and `.lmebayes_mer_optional_args()` are part of
the shared front door and are documented in
`inst/README_LMERB_GLMERB_FRONT_DOOR.md`, not here.

### 2.3 C++ functions specific to this model (`lmebayesCore/src`)

| File | Function | Exported? | Role |
|---|---|---|---|
| `src/twoBlockGibbs.cpp` | `two_block_rNormal_reg_v2_cpp_export()` | Yes | Full two-block Gibbs driver: `n` outer draws x `m_convergence` inner sweeps, alternating Block 1 / Block 2. Model-specific logic here: the two-block bookkeeping (below), not the per-group draw itself. |
| `src/twoBlockGibbs.cpp` | `block1_prior_list(...)` | No (internal) | Rebuilds the Block-1 prior list (per-group `mu` from `mu_all`) at each sweep -- this model's `gamma -> mu_all -> Block 1 prior` coupling. |
| `src/twoBlockGibbs.cpp` | `block2_prior_prep_v2(...)` / `struct Block2PriorV2` | No (internal) | Parses one Block-2 (`dNormal`, here) pfamily component once per outer draw. |
| `src/progress_utils.cpp` | `progress_bar(...)` / `progress_bar_finish(...)` | No (internal; via `glmbayes::progress::...`) | Console progress-bar rendering during the outer draw loop. |

The per-group conjugate-Normal draw itself (`block_rNormalReg_cpp_export()` ->
`rNormalRegBlocks()` -> `rNormalReg()`) and the `mu_all` builder
(`MuAllBuilder`/`build_mu_all()`) are shared dependencies -- see Section 3.

---

## 3. Shared dependencies (out of scope here; tracked separately)

This model's Gibbs engine calls into generic infrastructure that is not
specific to fixed-dispersion/known-vcov LMMs and should be evaluated for
porting on its own terms, separately from this document:

| Function | File | Why out of scope |
|---|---|---|
| `normalize_block()` (R) / `normalize_block_cpp()` (C++) | `R/simfunction_block_utils.R`, `src/block_utils.cpp` | Generic row-partition normalization; also the function `lmbBlock()`/`glmbBlock()` call directly. |
| `rNormalReg()` | `src/rNormalReg.cpp` | Generic conjugate Normal-regression draw; no LMM-specific (`gamma`/`tau2`/ICM) logic. |
| `rNormalRegBlocks()` | `src/rNormalRegBlocks.cpp` | Generic "loop `rNormalReg()` over row-blocks." |
| `block_rNormalReg_cpp_export()`, `normalize_prior_for_blocks_cpp()`, `prior_payload_from_blocks()` | `src/block_utils.cpp` | Generic per-block prior normalization/packing around `rNormalRegBlocks()`. |
| `.rNormalReg_cpp()` / `.block_rNormalReg_cpp()` / `.rNormalRegBlocks_cpp()` (R wrappers) | `R/rcpp_wrappers.R` | Thin `.Call()` bridges for the above. |
| `build_mu_all()` / `build_mu_all_r()` / `.lmerb_validate_design()` | `R/build_mu_all.R` | Maps a hyper design + coefficient vector to a group-level linear predictor; used by this model but generic beyond it. |
| `MuAllBuilder`, `two_block_build_mu_all()`, `two_block_build_mu_all_cpp_export()` | `src/two_block_block1.cpp`, `src/export_wrappers.cpp` | C++ backend for `build_mu_all()`. |
| `.two_block_build_mu_all_cpp()` (R wrapper) | `R/rcpp_wrappers.R` | Thin `.Call()` bridge for the above. |

No further porting-priority discussion for these is included here; that
decision spans other routes and callers (e.g. `lmbBlock()`/`glmbBlock()`) and
belongs in its own document.

---

## 4. Verification method / caveats

Every function/file pairing above was confirmed against the source in this
workspace (`c:\Rpackages\lmebayesCore`, `c:\Rpackages\lmebayes`) by direct
reading and `grep`, not inferred from documentation or memory; in particular:

- The C++ call chain (Section 1.2) was confirmed by reading
  `R/rcpp_wrappers.R` (positional `.Call()` targets) and
  `src/twoBlockGibbs.cpp` (the `v2` driver body and its direct calls to
  `block_rNormalReg_cpp_export()` for Block 1 and `rNormalReg()` for Block 2).
- The shared-dependency boundary (Section 3) was confirmed by reading
  `src/block_utils.cpp`, `src/rNormalRegBlocks.cpp`, `src/rNormalReg.cpp`,
  `src/two_block_block1.cpp`, and `R/build_mu_all.R`.
- The front-door split (Section 1 header, Section 2.1/2.2 pointers) was
  confirmed by reading `R/mixed_rmerb_helpers.R`, `R/model_setup.R`,
  `R/lmebayes_reg_route_table.R`, `R/rlmerb.R`, and `R/rglmerb.R`; see
  `inst/README_LMERB_GLMERB_FRONT_DOOR.md` for that half of the call graph.
- This document reflects the code as of the `sim_method` feature addition
  (see `NEWS.md` in both packages); if the two-block Gibbs C++ driver is
  later changed to call a different export (e.g. a `v3`), Sections 1.2, 2.3,
  and 3 will need to be re-verified against `R/rcpp_wrappers.R` and
  `src/twoBlockGibbs.cpp`.
