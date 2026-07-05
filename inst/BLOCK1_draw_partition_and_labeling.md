# Block 1 draw: row partition, prior wiring, and coefficient labeling

Maintainer-facing reference for how **conditionally independent block simulation**
partitions observations, assigns priors to blocks, and labels returned coefficients.
Use this when implementing or reviewing a **C++ batch loop** that replaces the
production R draw path (`.two_block_block1_draw_all_chains_block_outs()`).

Related docs:

- `inst/DESIGN_RGLM_BLOCKS.md` — block GLM design and return layout
- `inst/DESIGN_RNORMALREG_BLOCKS.md` — Gaussian block analogue
- `inst/PLAN_block1_cpp_migration.md` — migration plan and call tree

---

## Production path today (R draw loop)

Ex_16 / `rGLMM_sweep` with `R_engine`:

```
rGLMM_sweep
  └── .two_block_block1_all_chains()
        ├── .two_block_block1_prep_all_chains()     # mu_all, prior_list per chain
        └── .two_block_block1_draw_all_chains()
              ├── .two_block_block1_draw_all_chains_block_outs()   # for-loop over chains
              │     └── .two_block_block1_draw_block()
              │           └── block_rNormalGLM() / block_rNormalReg()   # R wrappers
              └── .two_block_block1_draw_reorder() per chain
```

Each chain draw ultimately reaches the same C++ backend as a direct export call:

```
block_rNormalGLM()  →  .block_rNormalGLM_cpp()
                    →  block_rNormalGLM_cpp_export()   # src/block_utils.cpp
                    →  rNormalGLMBlocks()               # src/rNormalGLMBlocks.cpp
                    →  rNormalGLM() per block
```

The R wrappers add **dimnames** after the `.Call`; reorder runs in R (or C++ helper)
after the draw.

---

## Step 1 — Build the partition (`normalize_block`)

**R:** `normalize_block(block, l2)` in `R/simfunction_block_utils.R`  
**C++:** `normalize_block_cpp()` / `split_factor_rows()` in `src/block_utils.cpp`

For the usual case — `design$groups` is a **factor** of length `l2`:

| Field | Meaning |
|-------|---------|
| `k` | Number of blocks (levels) |
| `ids[j]` | Block label = `levels(groups)[j]` |
| `rows[[j]]` | **1-based** row indices into full `y`, `Z`, `offset`, `weights` for block `j` |
| `l2_blocks[j]` | `length(rows[[j]])` |
| `starts[j]` | Cumulative start index in stacked layout (legacy; slicing uses `rows`) |

**Factor level order is preserved.** For a factor, `ids` follows `levels(block)`, not
sorted labels and not order of first appearance in the data. This matters for matching
`batch$group_levels` and for reorder (below).

---

## Step 2 — Assign prior mean column to block (`normalize_prior_for_blocks`)

**R:** `normalize_prior_for_blocks()` in `R/simfunction_block_utils.R`  
**C++:** `normalize_prior_for_blocks_cpp()` in `src/block_utils.cpp`

Block 1 passes `prior_list$mu` as an **`l1 × k`** matrix (`l1 = ncol(Z)`, `k` blocks).
When priors differ by block, column **`j` is assigned to block `j`** by **position**:

```r
for (j in seq_len(k)) {
  pl_j <- list(mu = mu[, j], ...)
}
```

**There is no lookup by `colnames(mu)` inside the block simulators.** Column `j` must
be the prior mean for `block_info$ids[j]`.

`build_mu_all()` produces columns in **`group_levels`** order (see `R/build_mu_all.R`).
In the standard `rGLMM_reg` setup, `group_levels` and `levels(design$groups)` are the same
vector, so `mu[, j]` aligns with `ids[j]`. If those orderings ever differ, prep must
permute `mu_all` columns to `block_info$ids` order **before** draw — reordering
coefficient rows after draw cannot fix a wrong prior assignment.

`prior_payload_from_blocks()` then stacks per-block priors into the payload consumed by
`rNormalGLMBlocks()` (`mu` matrix, `P_blocks`, `prior_by_block` flag).

---

## Step 3 — Subset rows per block (`rNormalGLMBlocks`)

**File:** `src/rNormalGLMBlocks.cpp` (Gaussian: `src/rNormalRegBlocks.cpp`)

For each block index `b` (0-based in C++, block `b+1` in R):

```cpp
IntegerVector rows = row_blocks[b];
y_b  = slice_numeric(y, rows);
x_b  = slice_matrix_rows(x, rows);
offset_b = slice_numeric(offset_full, rows);
wt_b     = slice_numeric(wt_full, rows);
mu_b     = mu_for_block(mu, b, prior_by_block);   // column b when prior_by_block
P_b      = P_for_block(P_blocks, b, prior_by_block);
```

- Slicing uses **1-based** R indices from `block_info$rows`.
- One **`rNormalGLM()`** / **`rNormalReg()`** call runs on the sub-problem only.
- Block Gibbs production uses **`n = 1`** per chain per sweep.

**Offset / weights:** `block_rNormalGLM()` defaults to `rep(0, l2)` and `rep(1, l2)`
when `offset` / `weights` are omitted. The Block 1 R draw wrapper does not pass
`design$offset` or `design$weights`. Any C++ orchestrator must match that contract
unless the R wrapper is changed deliberately.

---

## Step 4 — Stack block draws into `coefficients`

Inside `rNormalGLMBlocks()`:

| Output | Layout |
|--------|--------|
| `coefficients` | **`k × l1`** matrix |
| `coef.mode` | **`k × l1`**, same layout |
| Row `b` | Draw for block `b` (= `block_info$ids[b+1]`) |
| Col `c` | Coefficient `colnames(Z)[c]` |

The C++ export returns these matrices **without dimnames**. Row order matches
**block index order** from `normalize_block`, i.e. **`block_info$ids` order**.

---

## Step 5 — R wrapper labels (`block_rNormalGLM` / `block_rNormalReg`)

**File:** `R/simfunction_block.R`

After `.block_*_cpp()` returns:

```r
colnames(coef_draw) <- colnames(x)      # colnames(Z)
rownames(coef_draw) <- block_info$ids   # state / group labels
```

Same for `coef.mode`. The returned list includes `block_info` for downstream use.

---

## Step 6 — Reorder to `batch$group_levels` (Block 1 batch path)

**File:** `R/two_block_batch_gibbs.R` — `.two_block_block1_draw_reorder()`

Storage in `batch$b` uses row order **`batch$group_levels`**. Reorder permutes rows of
the draw matrix:

```r
ord <- match(group_levels, rownames(b_draw))
b_draw <- b_draw[ord, , drop = FALSE]
```

Requires **`rownames(b_draw)`** set to block ids (from Step 5). If rownames are
missing, R reorder is a no-op and rows may be misaligned with `group_levels` when
those orderings differ.

C++ helper: `two_block_reorder_b_to_group_levels()` (used in
`two_block_block1_one_chain_from_mu_P_impl()` in `src/two_block_block1.cpp`).

---

## R↔C++ carry-over checklist

C++ moves bare matrices and lists. R paths rely on **dimnames**, **names** on
vectors/lists, **factor attributes**, and **positional pairing** (column `j` ↔ block
`j` ↔ `ids[j]`). None of that is enforced inside the block simulator — each boundary
must copy or reconstruct metadata explicitly.

Use this checklist when moving Block 1 orchestration **R → C++ → R** (prep, draw,
reorder, batch assign). Items marked **silent** fail without error when wrong.

### Handoff A — Batch state → prep

| Must carry | Where set (R) | Risk if lost in C++ |
|------------|---------------|---------------------|
| `fixef[[k]]` colnames → `names(fixef_i)` | `.two_block_batch_fixef_chain()` | `build_mu_all` / `X_hyper` misalignment |
| `tau2[i, ]` colnames = `re_names` | batch init | ING `P` refresh uses wrong component index |
| `batch$b` dimnames `(group_levels, re_names, NULL)` | `.rGLMM_sweep_initialize()` | Unnamed 3-D array after slice assign |
| `group_levels` vs `levels(design$groups)` | `rGLMM_reg()` | Divergent order breaks mu, reorder, storage |

**Verify:** `names(fixef_i[[k]])` match `colnames(X_hyper[[k]])`; `names(tau2_i)` match
`re_names`.

### Handoff B — Prep output → draw input

| Must carry | Notes | Risk |
|------------|-------|------|
| `mu_all` rownames = `re_names`, colnames = `group_levels` | Set by `build_mu_all` / C++ `set_matrix_dimnames` | Diagnostics only unless prep permutes columns |
| `prior_list$mu[, j]` ↔ `block_info$ids[j]` | **By column index, not `colnames(mu)`** | **Wrong prior per state (silent)** |
| `prior_list$P`, `dispersion`, `ddef` | List elements round-trip | Wrong envelope / ING precision |
| Column order = `ids` order | Same set as `group_levels` may differ in order | Reorder after draw **cannot** fix wrong mu |

**Verify:** For each `j`, `mu_all[, j]` is the prior mean for state `block_info$ids[j]`
(e.g. compare to `build_mu_all` indexed by name, not just `colnames(mu)[j]`).

### Handoff C — `design` → C++ export

| Must carry | Notes | Risk |
|------------|-------|------|
| `design$groups` as **factor** SEXP | Levels attribute preserved | Re-coercion sorts labels; `ids` order changes |
| `colnames(Z)` | Not passed as separate arg | Coef columns unnamed after export |
| Default `offset` / `weights` | R wrapper: `rep(0,l2)`, `rep(1,l2)` | Reading `design$offset` changes likelihood **silent** |
| `family` + `link` strings, `f2`/`f3` | R resolves via `glmbfamfunc()` | Wrong sampler if strings or closures stale |

**Verify:** `identical(block_info$ids, levels(design$groups))` for factor groups.

### Handoff D — C++ export → labeled coefficients (**largest gap**)

`block_rNormalGLM_cpp_export()` returns **`k × l1` matrices without dimnames** plus
`block_info`. The R wrapper in `block_rNormalGLM()` attaches labels:

```r
colnames(coef_draw) <- colnames(x)
rownames(coef_draw) <- block_info$ids
```

| If skipped in C++ path | Effect |
|------------------------|--------|
| `rownames(coefficients)` | Reorder no-op or wrong (**silent** in R when NULL) |
| `colnames(coefficients)` | `batch$b` columns lose RE names |
| `block_info` on return | Reorder / diagnostics lose partition |
| `block_results` | `iters_mean` → `1` (**silent**) |
| S3 `class` | Not needed for Gibbs |

**Verify:** After draw, before reorder:
`identical(rownames(coefficients), block_info$ids)` and
`identical(colnames(coefficients), colnames(design$Z))`.

C++ reference: `set_block_draw_coefficient_dimnames()` in `two_block_block1.cpp`
(mirrors R wrapper).

### Handoff E — Draw output → reorder → `batch$b`

| Must carry | Notes | Risk |
|------------|-------|------|
| Row order before reorder | `block_info$ids` order | Coefficient row `j` = block `j` |
| `group_levels` target order | Separate arg to reorder | Permutation via `match(group_levels, rn)` |
| Row/col dimnames on reordered `b` | C++ sets in `ensure_batch_b_dimnames` | Slice into `batch$b` mislabels chains |

**Verify:** `identical(rownames(b_reordered), group_levels)`; `identical(colnames(b_reordered), re_names)`.

**Silent failure:** R `.two_block_block1_reorder_b_r()` and C++
`two_block_reorder_b_to_group_levels()` both **return unchanged** when rownames /
`block_ids` are NULL.

### Handoff F — `block_results` → `iters_ranef`

| Must carry | Notes | Risk |
|------------|-------|------|
| `block_out$block_results` | List length `k`, each with `iters` | Missing → `iters_mean = 1` (**silent**) |
| Mean over blocks | `.two_block_block1_iters_mean()` | Under/over-count envelope work |

**Verify:** `length(block_out$block_results) == block_info$k` after draw.

### Handoff G — Mixed piecewise C++ flags

When only some steps use C++ (`use_cpp_mu_all`, `use_cpp_prior_tau2`,
`use_cpp_reorder`, `use_cpp_iters`), numeric values may match while **names/order**
diverge. Production `rGLMM_sweep` defaults: prep/reorder/iters **R** (`FALSE` flags).

| Flag | Default (sweep path) | Carry-over concern |
|------|----------------------|-------------------|
| `use_cpp_mu_all` | `FALSE` | C++ must set same `mu_all` dimnames |
| `use_cpp_prior_tau2` | `FALSE` | ING `P` structure identical |
| `use_cpp_reorder` | `FALSE` | Same `match(group_levels, ids)` |
| `use_cpp_iters` | `FALSE` | Same `block_results` walk |

### Pre-flight checklist (run before trusting a C++ batch loop)

- [ ] **1.** `levels(design$groups)` and `batch$group_levels` agree (order per your contract)
- [ ] **2.** `prior_list$mu[, j]` is the mean for `block_info$ids[j]` (positional)
- [ ] **3.** `rownames(coefficients) == block_info$ids` after draw, before reorder
- [ ] **4.** `colnames(coefficients) == colnames(Z)` after draw
- [ ] **5.** `rownames(b)` after reorder equals `group_levels` in order
- [ ] **6.** `dimnames(batch$b)` is `(group_levels, re_names, NULL)` after assign
- [ ] **7.** `block_results` present when accumulating `iters_ranef`
- [ ] **8.** Offset `0`, weights `1` unless R wrapper contract changes
- [ ] **9.** `design$groups` passed as factor, not integer codes without levels
- [ ] **10.** `names(fixef_i[[k]])` and `names(tau2_i)` intact through prep

### Does not need to carry (common confusion)

- **R RNG state** — not shared; not a carry-over bug, but invalidates seed-matched draw tests
- **`colnames(mu)` inside block sim** — ignored; only column **index** matters
- **S3 `class` on block draw** — batch path does not use it
- **Full `family` object through C++** — only `family`/`link` strings + `f2`/`f3` at export

---

## Implications for a C++ version of the all-chains draw loop

The function under development mirrors
**`.two_block_block1_draw_all_chains_block_outs()`** — a loop over chains calling the
block simulator once per chain. Existing C++ pieces include:

| Piece | Location | Role |
|-------|----------|------|
| `two_block_block1_one_chain_from_mu_P_impl()` | `src/two_block_block1.cpp` | Prep + single draw + dimnames + reorder (reference for one chain) |
| `block_rNormalGLM_cpp_export()` | `src/block_utils.cpp` | Same backend as `.block_rNormalGLM_cpp()` |
| `two_block_block1_draw_block_impl()` | (experimental / reverted) | Direct export call + post-draw dimnames |

### Must match the validated R path

1. **Backend** — Call `block_rNormalGLM_cpp_export` / `block_rNormalReg_cpp_export`
   with the **same arguments** the R wrappers pass to `.block_*_cpp`, not a different
   code path or R callback to `block_rNormalGLM()`.

2. **Partition** — Pass `design$groups` as `block`; let `normalize_block_cpp` build
   `rows` and `ids`. Do not reimplement partition logic with different level ordering.

3. **`prior_list$mu` columns** — Column `j` must correspond to `block_info$ids[j]`.
   Prep (`build_mu_all` + `prior_with_tau2`) must agree with partition order before
   the C++ loop runs.

4. **Offset / weights** — Use `rep(0, l2)` and `rep(1, l2)` when the R wrapper omits
   them. Do not read `design$offset` / `design$weights` unless the R contract changes.

5. **Dimnames after export** — C++ export returns unnamed `k × l1` matrices. Mirror
   the R wrapper: `colnames <- colnames(Z)`, `rownames <- block_info$ids`, before
   reorder. See `set_block_draw_coefficient_dimnames()` in `two_block_block1.cpp`.

6. **Reorder** — After draw, permute rows to `batch$group_levels` (same as
   `.two_block_block1_draw_reorder()`). Reorder fixes **row storage order**, not
   wrong **`mu` column assignment**.

7. **GLM args** — Match R wrapper: `use_parallel = FALSE`, `verbose = FALSE`,
   `n_envopt = 1` for Block 1 Gibbs (`n = 1`).

### Do not require RNG parity between R and C++ loops

The envelope sampler does not use R's RNG. **Do not** validate a C++ batch loop by
matching draws to the R loop with `set.seed()`. Test wiring with **fixed** inputs
(made-up `mu_all`, `prior_list`, fake coefficient matrices) and compare
**distributions** on full runs (fixef means, state tables), not draw-by-draw equality.

### Suggested deterministic checks (no simulation RNG)

See **R↔C++ carry-over checklist** above for the full handoff map. Minimal tests:

| Check | What to verify |
|-------|----------------|
| Partition | `identical(block_info$ids, levels(design$groups))` and rows cover `1:l2` |
| Prior columns | For each `j`, `mu[, j]` equals `build_mu_all` value for `ids[j]` |
| Dimnames | `rownames(coefficients)` equals `block_info$ids` after post-export step |
| Reorder | `rownames(reordered_b)` equals `group_levels` in order |
| Slice size | Sum of `length(rows[[j]])` equals `l2`; each block draw uses only those rows |

Use **fabricated** `mu_all`, `prior_list`, and coefficient matrices — not two RNG draws.

### Anti-patterns observed in failed C++ batch experiments

- Calling `design_offset_wt(design, …)` instead of default zero/one vectors when R
  wrappers omit offset/weights.
- Skipping R-wrapper dimnames and relying on reorder alone (rownames may be missing
  or misaligned).
- Assuming `build_mu_all` column order always equals block index order without checking
  `group_levels` vs `block_info$ids`.
- Replacing the export backend with R callbacks in the hot loop (correct but too slow).
- Renaming or bypassing `*_cpp_export` entry points instead of calling the same
  symbols `.block_rNormalGLM_cpp` uses.

---

## Quick reference diagram

```
design$groups (factor, length l2)
    │
    ▼
normalize_block  ──►  block_info$ids[j]     block_info$rows[[j]]
    │                      │                        │
    │                      │                        └──► y[rows], Z[rows, ]  per block
    │                      │
prior_list$mu[, j] ────────┘  (by column index j, not name)

rNormalGLMBlocks  ──►  coefficients[k × l1]   (row b = block b, unnamed)

block_rNormalGLM (R)  ──►  rownames = ids,  colnames = colnames(Z)

draw_reorder  ──►  rows permuted to batch$group_levels  ──►  batch$b[, , chain]
```

---

## Files to read when changing this path

| Topic | File |
|-------|------|
| Partition (R) | `R/simfunction_block_utils.R` — `normalize_block()` |
| Prior split (R) | `R/simfunction_block_utils.R` — `normalize_prior_for_blocks()` |
| R draw wrapper | `R/simfunction_block.R` — `block_rNormalGLM()`, `block_rNormalReg()` |
| Export + payload | `src/block_utils.cpp` — `block_rNormalGLM_cpp_export()`, `prior_payload_from_blocks()` |
| Row slice loop | `src/rNormalGLMBlocks.cpp`, `src/rNormalRegBlocks.cpp` |
| R batch draw loop | `R/two_block_batch_gibbs.R` — `.two_block_block1_draw_all_chains_block_outs()` |
| C++ one-chain reference | `src/two_block_block1.cpp` — `two_block_block1_one_chain_from_mu_P_impl()` |
| Design spec | `inst/DESIGN_RGLM_BLOCKS.md` |
