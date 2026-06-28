# Block 1 all-chains step C: `batch$b[, , i] <- out$b`

Reference: committed R loop (`two_block_block1_one_chain_cpp`, step C only).

## What R does

```r
batch$b[, , i] <- out$b
```

- `batch$b` has `dim = c(J, p_re, n)` from `.two_block_batch_init`.
- `out$b` is `J × p_re` (rows = groups in `group_levels` order after reorder).
- Assignment is **positional**: row `g`, col `j` of `out$b` → `batch$b[g, j, i]`.
- R does **not** match by rownames/colnames on subassignment; slice dimnames stay those of the array.

## C++ replication (`batch_b_assign_slice`)

Flat index for R column-major `dim = c(J, p_re, n)`:

```
b[g, j, chain_i]  →  b_store[g + J * (j + p_re * (chain_i - 1))]
```

with `g ∈ [0, J-1]`, `j ∈ [0, p_re-1]` (0-based).

Same formula as v5 `store_b_chain` in `twoBlockGibbs.cpp`.

## Likely divergence sources (all-chains path)

1. **Wrong flat index** (swapped `J`/`p_re`, 0/1-based chain index, row-major assumption).
2. **Lost `dim` attribute** on SEXP → wrong `J`, `p_re` from `attr(b, "dim")`.
3. **`array(out$b, dim=...)` rebuild** after return (removed from R driver; prefer `batch$b <- out$b` with preserved `dim`).
4. **`out$b` orientation** not `J × p_re` (caught by dimension check if `J ≠ p_re`).
5. **Confusing rownames with position** — reorder must finish before step C so row `g` is group `group_levels[g]`.

## Verification (values only)

Run `data-raw/test_block1_b_assign_slice_cpp.R`: random `b_store`, random `b_draw`, every chain index; `identical` R subassignment vs `two_block_batch_b_assign_slice_cpp_export`.

Do **not** compare single draws between paths (RNG differs). Compare long-run colMeans after enabling C++ step C.

## Wiring

| Control | Default | Effect |
|---------|---------|--------|
| `options(glmbayesCore.use_cpp_b_slice = FALSE)` | — | Step C stays pure R |
| `use_cpp_b_slice = TRUE` | — | Per-call C++ slice |
| `use_cpp_block1_all_chains = TRUE` | `FALSE` | Uses C++ step C inside `all_chains_impl` (same `batch_b_assign_slice`) |

Default production path: R loop + **R step C** until slice test passes, then opt in C++ step C.
