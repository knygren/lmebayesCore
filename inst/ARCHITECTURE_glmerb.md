# glmerb / two-block Gibbs — internal architecture

Maintainer-facing map of the **pure-R sweep-outer** path (`rGLMM_sweep`).
This is the reference implementation for two-block GLMM Gibbs sampling in the glmbayes
ecosystem. The legacy C++ v5 driver (`two_block_rNormal_reg_v5`) was written to mirror
the same loop order but is a separate port.

Excluded from the source tarball via `.Rbuildignore` (same convention as
`inst/DESIGN_RGLM_BLOCKS.md`).

---

## High-level flow (glmerb, R engine)

```
glmerb (lmebayes)
  └── rglmerb
        └── [R_engine] .lmebayes_run_glmm_engine()
              └── REG_ROUTE_TABLE
                    ├── rGLMM_reg_known_vcov()      # all Block~2 dNormal; standard eigenvalue bounds
                    │     └── pilot (non-Gaussian always unless n_pilot = 0) + rGLMM_sweep
                    └── rGLMM_reg_estimated_vcov()  # ING Block~2; disp_lower eigenvalue bounds
                          └── pilot + rGLMM_sweep
                    └── rGLMM_sweep     # outer sweep loop (main and pilot stages)
                          └── two_block_batch_gibbs.R # Block 1 / Block 2 batch updates
```

Default **`glmerb`** uses **`R_engine`** → `.rglmerb_v6_rGLMM` → **`.lmebayes_run_glmm_engine()`**
→ **`REG_ROUTE_TABLE`** → **`rGLMM_reg_*`** (pure-R replicate-chain orchestration). Non-Gaussian models **always** run a pilot
stage unless **`n_pilot = 0L`**; the two routes differ in **eigenvalue-bound complexity**
(standard fixed-dispersion vs ING **`disp_lower`** conservatism), not in whether a
pilot runs. Set `.rglmerb_engine <- "cpp_engine"` in
`lmebayes/R/rglmerb.R` to use `rglmerb_v5` → `two_block_rNormal_reg_v5` (legacy C++).

---

## Outer sweep (one stage: pilot or main)

| File | Function | Role |
|------|----------|------|
| `glmbayesCore/R/rGLMM_sweep.R` | `rGLMM_sweep()` | Outer sweep loop |
| `glmbayesCore/R/two_block_batch_gibbs.R` | `.rGLMM_sweep_initialize`, `.rGLMM_sweep_save` | Batch state init / pack draws |
| `glmbayesCore/R/two_block_sweep_history.R` | `.two_block_build_sweep_history`, etc. | Per-sweep fixef colMean tables |

Loop order is **sweep-outer**: for each inner sweep `m`, all chains run Block 1, then all
chains run Block 2.

```r
for (m in seq_len(inner_sweeps)) {
  batch <- .two_block_block1_all_chains(...)
  batch <- .two_block_block2_all_chains(...)
  sweep_stats[[m]] <- .two_block_snapshot_fixef_stats(batch, re_names)
}
```

After all sweeps, `.rGLMM_sweep_save()` produces `fixef_draws`, `coefficients`
(b), `sweep_history`, etc.

---

## Block 2 call chain (one sweep)

```
rGLMM_sweep
  └── .two_block_block2_all_chains          # for i in 1:n; optional progress bar
        └── .two_block_block2_one_chain     # core Block 2 logic (per chain)
              ├── .two_block_align_b_to_xhyper   # reorder b to X_hyper rows
              └── rglmb(..., pfamily = pf)       # one call per RE component
                    └── simfunction → rNormal_reg / rIndepNormalGamma_reg / …
              └── .two_block_rglmb_iter_count    # ING only
```

**Substance:** `.two_block_block2_one_chain` (~50 lines): align `b`, call `rglmb`, write
`batch$fixef[i, ]` from **`coefficients[1, ]`** (and `tau2` / `iters` if ING).

**`.two_block_block2_all_chains`:** loop over chains + optional progress bar only.

| File | Function | Role |
|------|----------|------|
| `glmbayesCore/R/rGLMM_sweep.R` | `rGLMM_sweep` | Outer `m` loop; calls Block 1 then Block 2 |
| `glmbayesCore/R/two_block_batch_gibbs.R` | `.two_block_block2_all_chains` | Loop `i = 1..n` |
| `glmbayesCore/R/two_block_batch_gibbs.R` | `.two_block_block2_one_chain` | Per-chain Block 2 update |
| `glmbayesCore/R/two_block_batch_gibbs.R` | `.two_block_align_b_to_xhyper` | Map `b` (`group_levels` order) → `y` (`X_hyper` row order) |
| `glmbayesCore/R/rglmb.R` | `rglmb` | Block 2 Gaussian regression on `pfamily` |
| `glmbayesCore/R/simfunction.R` | `rNormal_reg`, `rIndepNormalGamma_reg`, … | Draw / `coef.mode` via `pfamily$simfun` |
| `glmbayesCore/R/two_block_batch_gibbs.R` | `.two_block_rglmb_iter_count` | ING envelope iteration count |

**Align is critical:** `y_k` passed to `rglmb` must follow `rownames(X_hyper[[k]])`, with
`b` stored in `group_levels` order after Block 1. When rownames match `group_levels`
(e.g. book-banning), align is identity; other datasets expose mis-order bugs if align is
skipped.

---

## Block 1 call chain (one sweep)

Block 1 is split into **prep** then **draw** (optional parallel over chains).

```
rGLMM_sweep
  └── .two_block_block1_all_chains
        ├── .two_block_block1_prep_all_chains
        │     └── .two_block_block1_prep_one_chain
        │           ├── .two_block_batch_fixef_chain
        │           ├── build_mu_all
        │           └── .two_block_block1_prior_with_tau2
        └── .two_block_block1_draw_all_chains
              └── .two_block_block1_draw_one_chain
                    └── block_rNormalGLM / block_rNormalReg
```

**Substance:** `.two_block_block1_draw_one_chain` → `block_rNormalGLM` (GLMM) or
`block_rNormalReg` (Gaussian), then reorder `b` to `group_levels`.

| File | Function | Role |
|------|----------|------|
| `glmbayesCore/R/two_block_batch_gibbs.R` | `.two_block_block1_all_chains` | Prep + draw all chains |
| `glmbayesCore/R/two_block_batch_gibbs.R` | `.two_block_block1_prep_one_chain` | fixef → `mu_all` → Block 1 `prior_list` |
| `glmbayesCore/R/build_mu_all.R` | `build_mu_all` | Hyperparameter mean matrix for Block 1 (C++ default) |
| `glmbayesCore/R/two_block_batch_gibbs.R` | `.two_block_block1_prior_with_tau2` | ING τ² → Block 1 `P` refresh (C++ default) |
| `glmbayesCore/R/two_block_batch_gibbs.R` | `.two_block_block1_draw_one_chain` | One Block 1 draw |
| `glmbayesCore/R/simfunction_block.R` | `block_rNormalGLM`, `block_rNormalReg` | Grouped RE draw (C++ under `.block_rNormalGLM_cpp`) |

---

## Layer summary

| Layer | Role |
|--------|------|
| `rGLMM_sweep` | Outer `m` loop |
| `*_all_chains` | Loop over `n` replicate chains + progress bars |
| `*_one_chain` | Actual Gibbs update for one chain |
| `build_mu_all`, align, prior helpers | Setup for Block 1 / Block 2 |
| `block_rNormalGLM`, `rglmb` | Where randomness and posterior modes happen |

Many function names exist for batch bookkeeping, pilot staging, and optional parallel prep;
**one sweep of sampling** is: for each chain, Block 1 = `block_rNormalGLM`, Block 2 =
`rglmb` per RE column.

---

## Staging / pilot / wrappers (not the sampler core)

| File | Role |
|------|------|
| `lmebayes/R/glmerb.R` | User API; `model_setup`, priors, calls `rglmerb` |
| `lmebayes/R/rglmerb.R` | Engine switch: `cpp_engine` vs `R_engine` |
| `lmebayes/R/rglmerb_v5.R` | Pilot/main planning; calls `run_short_chains_v5` (C++) |
| `glmbayesCore/R/rGLMM_reg.R` | Pilot/main orchestration via `rGLMM_sweep` (R engine); two route exports |
| `glmbayesCore/R/two_block_pilot_cost.R` | `n_pilot`, `m_convergence_pilot` planning |
| `glmbayesCore/R/two_block_glmm_pilot_helpers.R` | `fixef.init`, pilot UB, staged output names |
| `glmbayesCore/R/two_block_sweep_history.R` | Per-sweep fixef colMeans tables |

---

## Legacy C++ full sweep (optional / parity work)

| File | Role |
|------|------|
| `glmbayesCore/src/twoBlockGibbs.cpp` | `two_block_rNormal_reg_v5_cpp_export` — sweep-outer in C++ |
| `glmbayesCore/R/two_block_rNormal_reg_v5.R` | R wrapper → v5 C++ |
| `lmebayes/R/run_short_chains_v5.R` | Thin wrapper for `rglmerb_v5` |

Intended to mirror the R batch driver. Block 2 uses **`rglmb()$coefficients`** (posterior
draw), aligned with v5 `rNormalReg` / ING envelope output — not `coef.mode`.

Helpers already in C++ include `two_block_align_b_col_to_x_rows` and
`two_block_block2_rglmb_gamma` (align via R `.two_block_align_b_to_xhyper` + `rglmb`).

---

## Planned Block 1 C++ migration (random-effects `b`, all-chains)

**Block 1** updates random effects **`b`** (prep: `build_mu_all` + prior; draw:
`block_rNormalGLM` for Poisson/binomial GLMM, or `block_rNormalReg` on the Gaussian v6 path).
There is **`use_cpp_block1`** on `rGLMM_sweep` (default `TRUE`). Set
`use_cpp_block1 = FALSE` for the R prep/draw loop oracle with piecewise C++ helpers.

**Full plan:** [`inst/PLAN_block1_cpp_migration.md`](PLAN_block1_cpp_migration.md)

Summary:

| Phase | Deliverable | Eliminates |
|-------|-------------|------------|
| 0 | Fix ING `P` refresh in `.two_block_block1_prior_with_tau2` + C++ export | — ✅ |
| 1 | `two_block_block1_all_chains_cpp_export` (reuse `MuAllBuilder`, v5 ING refresh) | R prep/draw chain loops |
| 2 | Fuse prep + draw inside C++ chain loop | Intermediate `prior_lists` SEXP |
| 3 | Native fam dispatch (no R `f2`/`f3`) | **`glmbfamfunc`** per chain |
| 4 | C++ posterior mode in `rNormalGLM` | R **`optim()`** per group |
| 5 | Optional: conjugate Gaussian Block 1 on v6 | **`lm.fit`** in `rNormalReg` |

Legacy v5 already runs Block 1 inside one C++ sweep (`block_rNormalGLM_cpp_export` with
`f2`/`f3` passed once). v6 migration **extracts that body** into an all-chains export wired
from `.two_block_block1_all_chains`.

Fix Phase 0 (ING τ² ↔ `P`) before trusting parity on Ex_22 binomial ING.

---

## Planned Block 2 C++ migration (hyperparameters `γ`, all-chains, zero R callbacks)

**Status:** Phase 0 partial — `two_block_block2_one_chain_cpp` has native align but still
calls **`rglmb()`** per RE component; R **`for (i in 1:n)`** over chains remains.

**Full plan:** [`inst/PLAN_block2_cpp_migration.md`](PLAN_block2_cpp_migration.md)

Summary:

| Phase | Deliverable | Eliminates |
|-------|-------------|------------|
| 0 | One-chain parity gate | — |
| 1 | `two_block_block2_all_chains_cpp_export` | R chain loop + per-chain `.Call` |
| 2 | Direct `rNormalReg` / `rIndepNormalGammaReg` C++ (reuse `block2_prior_prep_v2`) | **`rglmb()`**, **`stats::gaussian()`** |
| 3 | `two_block_block2_dNormal_draw` (Armadillo conjugate) | **`lm.fit`** in dNormal Block 2 |
| 4 | ING draw path without **`optim()`** | R optimizer in ING Block 2 |
| 5 | Docs, static grep CI, remove A/B flags | — |

Legacy v5 already calls simfuncs directly (no `rglmb`) but skips align and still hits
`lm.fit` / `optim`. The v6 native path must **align like R** and reach **zero R callbacks**
on the Block 2 hot path.

Test align first (R vs C++ `y_k` exact), then fixef colMeans at large `n` with
`inner_sweeps = 1`.

---

## What one sweep stores

| Quantity | After 1 full sweep (`inner_sweeps = 1`) |
|----------|----------------------------------------|
| `fixef_draws` | Block 2 γ, one row per chain (after Block 2) |
| `coefficients` | Block 1 b from that sweep’s Block 1 |
| `sweep_history` sweep = 1 | Fixef chain colMeans (matches `colMeans(fixef_draws[[k]])` when `inner_sweeps = 1`) |

Pilot **`fixef.init`** = colMeans over **`n_pilot`** chains after **`m_convergence_pilot`**
sweeps (e.g. 19), not after one sweep.

For parity checks between R and C++ kernels: one `set.seed()` before each run, large
`n`, compare chain colMeans (fixef and b) — statistical agreement, not bit-exact RNG
matching.

---

## Related docs

- `inst/BLOCK_GIBBS_ERGODICITY.md` — ergodicity / full-rank conditions for block Gibbs
- `inst/DESIGN_RGLM_BLOCKS.md` — conditionally independent GLM block sampling
- `README.md` — package ecosystem and `rglmb()` orchestration
