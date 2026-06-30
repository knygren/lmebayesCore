# Plan: Block 1 all-chains C++ migration (random-effects updates)

Maintainer-facing plan to move **Block 1** — the **random-effects `b` draw** — off the
R batch loop used by `rGLMM_sweep` / non-Gaussian `glmerb` (`R_engine`).

**Terminology (two-block Gibbs):**

| Block | Updates | Typical call | Outer GLMM family? |
|-------|---------|--------------|-------------------|
| **Block 1** | Random effects **`b`** | `block_rNormalGLM` / `block_rNormalReg` | **Yes** (Poisson, binomial, …) for GLMM path |
| **Block 2** | Hyperparameters **`γ`** (+ **`τ²`** for ING) | `rglmb(..., gaussian())` | **No** — always Gaussian hyper-regression on `b` |

See `inst/ARCHITECTURE_glmerb.md` for the full call tree. Block 2 migration:
`inst/PLAN_block2_cpp_migration.md`.

**Out of scope:** Gaussian glmerb with all `dNormal` RE → `two_block_rNormal_reg_v2`
(full C++ sweep already). This plan targets the **v6 batch driver** Block 1 path.

---

## Goal

After migration, one inner sweep’s Block 1 phase is:

```
rGLMM_sweep
  └── .Call(two_block_block1_all_chains_cpp_export, ...)   # single entry
        └── C++ for (i) chains
              ├── MuAllBuilder.build(fixef_i)               # = build_mu_all
              ├── refresh ING rows of P from tau2_i         # = prior_with_tau2
              ├── block_rNormalGLM_cpp_export(1, ...)       # GLMM (Poisson/binomial/…)
              │     or block_rNormalReg_cpp_export(1, ...)  # Gaussian + ING on v6 path
              ├── reorder b to group_levels
              └── accumulate iters_ranef
```

**Non-goals:**

- Block 2 C++ migration (separate plan).
- Bit-exact RNG parity with R reference (statistical parity at large `n` only).
- Replacing legacy v5 driver — reuse its Block 1 logic where possible.

---

## Current state (R batch driver)

```
.two_block_block1_all_chains
  ├── .two_block_block1_prep_all_chains     # lapply / mclapply over chains
  │     └── .two_block_block1_prep_one_chain
  │           ├── .two_block_batch_fixef_chain
  │           ├── build_mu_all              # C++ default (use_cpp_mu_all)
  │           └── .two_block_block1_prior_with_tau2  # C++ default (use_cpp_prior_tau2)
  └── .two_block_block1_draw_all_chains     # lapply / mclapply over chains
        └── .two_block_block1_draw_one_chain
              └── block_rNormalGLM / block_rNormalReg
                    └── glmbfamfunc → .block_rNormalGLM_cpp → rNormalGLMBlocks → rNormalGLM
```

| Layer | Location | Already C++? | R callbacks? |
|-------|----------|--------------|--------------|
| Chain loop | R `lapply` / `mclapply` in prep + draw | No | R loop every sweep |
| `build_mu_all` | `R/build_mu_all.R` | **Yes** — `two_block_build_mu_all_cpp_export` | Optional R via `use_cpp_mu_all = FALSE` |
| ING `P` refresh | `.two_block_block1_prior_with_tau2` | **Yes** — `two_block_block1_prior_with_tau2_cpp_export` (v5 `any_ing`) | Optional R via `use_cpp_prior_tau2 = FALSE` |
| GLMM draw | `block_rNormalGLM()` | **Partial** — `block_rNormalGLM_cpp_export` | **`glmbfamfunc`** per chain; **`optim()`** per group inside `rNormalGLM.cpp` |
| Gaussian draw | `block_rNormalReg()` | **Partial** — `block_rNormalReg_cpp_export` | **`lm.fit`** / mode path in `rNormalReg.cpp` |
| `b` reorder | `match(group_levels, rownames)` | **Yes** — `two_block_reorder_b_to_group_levels_cpp_export` | Optional R via `use_cpp_reorder = FALSE` |
| `iters_ranef` mean | `.two_block_block1_iters_mean` | **Yes** — `two_block_block1_iters_mean_cpp_export` | Optional R via `use_cpp_iters = FALSE` |

### Incremental v6 exports (piece-by-piece)

| # | Piece | v5 C++ | v6 export | Done? | Notes |
|---|-------|--------|-----------|-------|-------|
| 1 | fixef extract | — | R only | No | `.two_block_batch_fixef_chain` |
| 2 | `build_mu_all` | `MuAllBuilder` | `two_block_build_mu_all_cpp_export` | **Yes** | `use_cpp_mu_all = TRUE` default |
| 3 | `prior_with_tau2` | inline ING loop | `two_block_block1_prior_with_tau2_cpp_export` | **Yes** | R `ddef` bug fixed; v5 `any_ing` semantics |
| 4 | Block 1 draw | `block_rNormalGLM_cpp_export` | via R wrappers | Partial | `glmbfamfunc` per chain remains |
| 5 | reorder `b` | v5 `block_info$ids` | `two_block_reorder_b_to_group_levels_cpp_export` | **Yes** | v6 `group_levels` semantics |
| 6 | `iters_mean` | inline in v5 | `two_block_block1_iters_mean_cpp_export` | **Yes** | `use_cpp_iters = TRUE` default |
| 7 | one-chain composite | prep + draw in R | `two_block_block1_one_chain_cpp_export` | **Yes** | `two_block_block1_one_chain_cpp()` |
| 8 | all-chains composite | R double loop | `two_block_block1_all_chains_cpp_export` | **Yes** | `use_cpp_block1` on `rGLMM_sweep` (default `TRUE`) |

There is **`use_cpp_block1`** on `rGLMM_sweep` (default `TRUE`; set `FALSE` for R prep/draw loops).

Legacy v5 (`twoBlockGibbs.cpp` inner sweep) already runs Block 1 in C++:

```cpp
mu_all = mu_builder.build(fixef);
pl1 = block1_prior_list(...);   // + ING P refresh when any_ing
block_i = block_rNormalGLM_cpp_export(1, y, x, block, pl1, ..., f2, f3, family, link, ...);
// or block_rNormalReg_cpp_export when is_gaussian
```

The v6 refactor **lifts this per-chain body** into an all-chains export wired from
`.two_block_block1_all_chains`.

---

## Correctness gate (before / alongside migration)

### ING prior precision coupling

Poisson/binomial glmerb uses `.lmebayes_block1_prior_list(..., dispersion_ranef = NULL)`
→ `ddef = TRUE`. R `.two_block_block1_prior_with_tau2` **used to** return early and skip
ING refresh; fixed to mirror v5 `any_ing` (refresh ING diagonal entries regardless of
`ddef`; still forward `ddef` for the envelope sampler dispersion convention).

**Fix in R first** (mirror v5): refresh ING diagonal entries regardless of `ddef`;
still forward `ddef` for the envelope sampler dispersion convention.

**Test:** `lmebayes/demo/Ex_22_glmerb_book_banning_ING.R`; extend
`lmebayes/data-raw/test_ing_sampling.R` for Block 1 τ² ↔ P on binomial path.

---

## R callback inventory (Block 1 hot path)

| Callback | Where | Needed? | Remove how |
|----------|-------|---------|------------|
| R `lapply` / `mclapply` over chains | prep + draw | No | Phase 1: all-chains export |
| `build_mu_all()` | prep | No | Done — `two_block_build_mu_all_cpp_export` |
| `.two_block_block1_prior_with_tau2` | prep | No | Done — `two_block_block1_prior_with_tau2_cpp_export` |
| `block_rNormalGLM()` R wrapper | draw | No | Phase 1: direct `.Call(block_rNormalGLM_cpp_export)` |
| `glmbfamfunc(family)` | inside `block_rNormalGLM()` | No | Phase 1: once per sweep; Phase 3: native `fam` |
| `f2` / `f3` R `Function` | `rNormalGLM` envelope | No | Phase 3: C++ fam dispatch |
| `optim()` | `rNormalGLM.cpp` posterior mode | Debatable | Phase 4: C++ mode finder |
| `block_rNormalReg()` R wrapper | Gaussian Block 1 on v6 | No | Phase 1: direct export |
| `lm.fit` | `rNormalReg.cpp` | For stored `b` mode | Phase 5 (Gaussian Block 1 only): conjugate draw |

---

## Target API

### New export (Phase 1)

```r
two_block_block1_all_chains_cpp_export(
  batch_b,              # J x p_re x n
  batch_fixef,          # list of n x q_k matrices
  batch_tau2,           # n x p_re
  batch_iters_ranef,    # length n
  y, Z, groups, offset, wt,
  x_hyper,              # for MuAllBuilder
  group_levels,
  re_names,
  ptypes,               # for ING P refresh
  block1_prior,         # base P, dispersion, ddef
  f2, f3,               # from one glmbfamfunc(family) call in R
  family, link,
  is_gaussian,
  Gridtype, n_envopt,
  use_parallel, use_opencl,
  progbar, progbar_prefix
)
# → list(b = ..., iters_ranef = ...)  same layout as current batch
```

R wrapper:

```r
.two_block_block1_all_chains(..., use_cpp_block1 = TRUE)
  → single .Call(two_block_block1_all_chains_cpp_export, ...)
```

Keep `.two_block_block1_all_chains` R path as **reference oracle** for parity tests.

---

## Phased work

### Phase 0 — ING `P` refresh fix (R) ✅

**Work:** Fix `.two_block_block1_prior_with_tau2` to match v5 `any_ing` logic; add
`two_block_block1_prior_with_tau2_cpp_export` with `use_cpp_prior_tau2` flag.

**Acceptance:** Ex_22 smoke; τ² from Block 2 changes Block 1 `P` on next sweep;
`data-raw/test_block1_prior_with_tau2_cpp.R` passes.

### Phase 1 — All-chains C++ loop (eliminate R chain loop)

**Work:**

1. Extract `two_block_block1_one_chain_impl(...)` from v5 inner loop (or port from
   `two_block_block1_prep_one_chain` + `two_block_block1_draw_one_chain` semantics).
2. Reuse **`MuAllBuilder`**, **`block1_prior_list`**, v5 ING **`P`** refresh.
3. Call **`block_rNormalGLM_cpp_export(1, ...)`** or **`block_rNormalReg_cpp_export(1, ...)`**
   per chain; pass **`f2`/`f3` once** from R (one `glmbfamfunc(family)` per Block 1 batch).
4. Port **`b` reorder** to `group_levels` (v5 reads `block_info$ids`; v6 uses explicit
   `match(group_levels, rownames)` — keep v6 semantics).
5. Port **`.two_block_block1_iters_mean`** logic for `iters_ranef`.
6. Add **`two_block_block1_all_chains_cpp_export`**; wire `use_cpp_block1` in
   `.two_block_block1_all_chains`.
7. Progress bar in C++ (reuse `glmbayes::progress::progress_bar`).

**R callbacks remaining:** `optim`, `f2`/`f3` inside `rNormalGLM` — **must not stop here**
if goal is zero callbacks.

**Acceptance:**

- dNormal Poisson/binomial: colMeans of `b` vs R reference at large `n`, `inner_sweeps = 1`.
- One `.Call` per Block 1 batch instead of `2n` (prep + draw) or `n` chain calls.
- Ex_16 / Ex_22 smoke.

### Phase 2 — Fuse prep + draw (no intermediate SEXP lists)

**Work:** Drop separate `prior_lists` / `mu_all` lists built for all chains in R; compute
`mu_all` and `pl1` inside the chain loop (as v5 does).

**Acceptance:** Same as Phase 1; lower allocation / SEXP traffic.

### Phase 3 — Native famfuncs (drop R `f2`/`f3`)

**Work:** Dispatch Poisson/binomial/Gamma links inside C++ (`glmbayes::fam` or equivalent);
stop passing R `Function` into `rNormalGLM`.

**Acceptance:** Grep Block 1 path: no `glmbfamfunc`, no `Function f2` in hot loop.

**Hurdle:** `DESIGN_RGLM_BLOCKS.md` notes fam dispatch is partially ported; link coverage
must match `block_rNormalGLM()` validation (poisson/log, binomial/logit|probit|cloglog, …).

### Phase 4 — Posterior mode without R `optim()` (GLMM Block 1)

**Work:** Replace `optim()` in `rNormalGLM.cpp` with C++ Newton/BFGS on the penalized
log-posterior (same as long-term goal for standalone `rglmb` GLM path).

**Acceptance:** Block 1 colMeans parity at large `n`; grep: no `Rcpp::Function optfun`.

**Hurdle:** Largest engineering item; envelope geometry depends on accurate mode + gradient.
Gridtype / link-specific edge cases (cloglog, probit).

### Phase 5 — Gaussian Block 1 on v6 path (optional)

When `is_gaussian` and model uses `rGLMM_sweep` (e.g. ING `lmerb`-style
routing), Block 1 calls `block_rNormalReg` instead of `block_rNormalGLM`.

**Work:** Conjugate / native mode in `rNormalReg` without `lm.fit` (similar to Block 2
Phase 3 for dNormal).

**Note:** All-`dNormal` Gaussian glmerb already uses `two_block_rNormal_reg_v2` — this
phase only matters for mixed or ING models on the v6 driver.

### Phase 6 — Docs, CI, flag cleanup

- Update `inst/ARCHITECTURE_glmbayes.md`.
- `data-raw/test_block1_all_chains_cpp.R` parity script.
- Optional: static grep CI for Block 1 R callbacks when `use_cpp_block1 = TRUE`.

---

## Major hurdles

### 1. ING τ² ↔ Block 1 `P` (correctness)

Block 2 draws τ²; Block 1 must see `P[k,k] = 1/tau2_k` on the next sweep. R v6 path
currently broken for `ddef = TRUE` (typical Poisson/binomial). **Fix before trusting any
C++ port.**

### 2. Two Block 1 simfuncs (routing)

| Condition | Draw function | Likelihood in draw |
|-----------|---------------|-------------------|
| `family != gaussian` (GLMM) | `block_rNormalGLM_cpp_export` | Poisson, binomial, … |
| `family == gaussian` (v6 only) | `block_rNormalReg_cpp_export` | Gaussian |

The all-chains export must branch on `is_gaussian` like `.two_block_block1_draw_one_chain`.

### 3. `rNormalGLM` still calls R `optim()` per group

Each chain × each grouping level runs envelope sampling after an **`optim()`** mode find.
This dominates cost for large `J`. Phase 4 is hard; Phase 1 still leaves this callback.

### 4. `f2` / `f3` R functions in the envelope sampler

Even with `block_rNormalGLM_cpp_export`, every `rNormalGLM` call uses R closures for
the GLM family. Phase 3 removes this; requires native link/family tables.

### 5. Group id ordering

Block partition returns coefficients with `block_info$ids` row order. v6 R path reorders
to **`batch$group_levels`**. C++ port must match v6 (not raw v5 `b_i` assignment).

### 6. `iters_ranef` bookkeeping

`.two_block_block1_iters_mean` averages envelope candidate counts across groups from
`block_results`. Easy to get wrong when collapsing to one chain update.

### 7. Parallelism

R prep/draw supports `mclapply` (not on Windows). v6 draw sets `use_parallel = FALSE`.
C++ `rNormalGLMBlocks` can parallelize over groups but batch driver disables it — decide
policy before enabling OpenMP inside all-chains export.

### 8. Parity standard

Match **chain colMeans** of `b` (and `iters_ranef` totals) at large `n`, not bit-exact
draws. Same convention as Block 2 plan.

### 9. Reuse vs duplication with v5

`twoBlockGibbs.cpp` already contains Block 1 inner logic embedded in full sweeps.
Prefer **extracting shared helpers** (`two_block_block1_one_chain_impl`) over a third
copy of ING refresh / reorder rules.

---

## Suggested order

1. **Phase 0** — ING `P` fix in R (small, blocks correct ING glmerb).
2. **Phase 1** — all-chains export + `use_cpp_block1` (structural win; reuses existing
   `block_rNormalGLM_cpp_export`).
3. **Phase 2** — fuse prep/draw in C++.
4. **Phase 3–4** — native fam + native optim (performance / zero callbacks).
5. **Phase 5** — only if v6 Gaussian Block 1 path needs it.

---

## Validation checklist

- [ ] Align / group order: `group_levels` vs `block_info$ids` on book-banning design.
- [ ] Poisson/binomial dNormal: R vs `use_cpp_block1` colMeans(`b`) at large `n`.
- [ ] Binomial ING (Ex_22): runs; τ² tracks Block 2; Block 1 `P` refresh each sweep.
- [ ] `data-raw/test_block_rNormalGLM_cpp.R` still passes (export unchanged).
- [ ] Profiling: one `.Call` per Block 1 batch vs current `O(n)` R loops.

---

## Related files

| Area | Path |
|------|------|
| R batch driver | `R/two_block_batch_gibbs.R` |
| R sweep | `R/rGLMM_sweep.R` |
| `build_mu_all` | `R/build_mu_all.R` |
| GLMM block draw | `R/simfunction_block.R`, `src/block_utils.cpp`, `src/rNormalGLMBlocks.cpp`, `src/rNormalGLM.cpp` |
| Gaussian block draw | `src/rNormalReg.cpp` |
| v5 reference loop | `src/twoBlockGibbs.cpp` (`MuAllBuilder`, Block 1 ~1640–1685) |
| Architecture | `inst/ARCHITECTURE_glmerb.md` |
| Block 2 plan | `inst/PLAN_block2_cpp_migration.md` |
