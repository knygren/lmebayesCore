# `lmbBlock()` — dependency assessment on `lmebayesCore`

Full call-graph trace of every function `lmbBlock()` (and its documented
alias `print.blmb`) actually invokes, isolating which of those are
**`lmebayesCore`** functions (R or C++) vs. functions from other packages.
Requested as a targeted companion to
[R_EXPORT_REACHABILITY.md](R_EXPORT_REACHABILITY.md), which already flags the
one-line summary this doc expands on (see its "lmebayes — additional entry
points" table).

**Scope note:** `lmbBlock()` itself is **not** defined in `lmebayesCore` — it
lives in the downstream package **`lmebayes`**
(`lmebayes/R/lmbBlock.R`). This doc traces *into* `lmebayesCore` from there.

**Last reviewed:** 2026-07-20 (against `lmebayes` `R/lmbBlock.R`,
`R/block_core_pfamily.R`; `lmebayesCore` `R/simfunction_block_utils.R`).

---

## Headline finding

`lmbBlock()` has exactly **one** runtime dependency on `lmebayesCore`:

```
lmebayesCore::normalize_block()
```

That's it. No other `lmebayesCore` R function, no `lmebayesCore` S3 method,
and **no `lmebayesCore` C++/compiled code** (`src/*.cpp`, `.Call()`) is
reached from `lmbBlock()`, directly or transitively. `normalize_block()`
itself is a self-contained, pure-R leaf function — it makes no further calls
into `lmebayesCore`, `glmbayesCore`, or any compiled code.

This is because `lmbBlock()` is a **`glmbayes::lmb()`-per-block** driver (SAS
`BY`-style independent fits), not a Gibbs/mixed-model sampler — the actual
posterior draws for each block happen entirely inside `glmbayes::lmb()`,
which is a fully self-contained package with its own compiled engine and
**zero** dependency on `lmebayesCore` (or `glmbayesCore`) at either the R or
C++ level. `lmbBlock()` only needs `lmebayesCore` to turn the user's `block`
argument into row-index groups before looping over `glmbayes::lmb()`.

---

## Call graph

```
lmbBlock()                                              [lmebayes/R/lmbBlock.R]
├── .blmb_formula_block_meta()                          [lmebayes, internal]
│   ├── stats::model.frame(), stats::model.matrix()     [base R]
│   ├── .blmb_resolve_block()                           [lmebayes, internal — pure R]
│   └── lmebayesCore::normalize_block()  ★ ONLY lmebayesCore CALL
│         └── (leaf: as.integer/is.list/lapply/unlist/duplicated/sort/
│              is.factor/nlevels/split/levels/vapply/cumsum/as.vector/factor
│              — all base R; no further lmebayesCore, glmbayesCore, or C++ calls)
│
├── .mrglmb_normalize_pfamily_lists()                   [from glmbayesCore, NOT lmebayesCore —
├── .validate_pfamily_for_rlmb()                          see "Look-alike dependencies" below]
│
├── (loop over k blocks)
│   └── do.call(lmb, ...)                               [glmbayes::lmb() — fully self-contained;
│                                                          zero lmebayesCore/glmbayesCore dependency]
│
├── .blmb_rows_to_data_subset()                          [lmebayes, internal — pure R]
├── .blmb_lmb_display_call()                             [lmebayes, internal — pure R]
└── .blmb_assemble()                                     [lmebayes, internal — pure R]

print.blmb(x, ...)                                       [alias documented on ?lmbBlock]
├── .blmb_coef_means_matrix()                            [lmebayes, internal — pure R]
└── .blmb_dic_table()                                    [lmebayes, internal — pure R]
```

No node in this graph other than `normalize_block()` resolves to
`lmebayesCore`. Verified by grepping `lmebayes/R/lmbBlock.R` for every
`lmebayesCore::` occurrence (two hits total — see next section) and reading
every function `lmbBlock()` transitively calls.

---

## The one dependency, in detail

### `lmebayesCore::normalize_block(block, l2)`

- **File:** [`R/simfunction_block_utils.R`](../R/simfunction_block_utils.R)
  (exported; `@export`, documented under `?block_simfuncs`).
- **Called from `lmbBlock()` via:** `.blmb_formula_block_meta()`, line
  `block_info <- lmebayesCore::normalize_block(block_vec, l2)`
  (`lmebayes/R/lmbBlock.R:180`).
- **Purpose:** normalizes the user-supplied `block` argument (factor,
  integer vector of length `l2`, `l2_blocks` counts, or a list of row-index
  vectors) into a canonical `list(k, ids, l2_blocks, starts, rows)` structure
  — `k` blocks, character `ids`, and `rows` (one integer index vector per
  block). `lmbBlock()` uses only `block_info$k` and `block_info$rows` (via
  `.blmb_rows_to_data_subset()`) from the returned list.
- **Implementation:** pure base R (`is.list`, `is.factor`, `split`,
  `vapply`, `cumsum`, `factor`, …). No `.Call()`, no Rcpp, no reference to
  any other `lmebayesCore` function, `glmbayesCore`, or compiled code.
  Confirmed by reading the full function body — it is a self-contained leaf.
- **Side effects:** none (pure function; only possible outcome besides a
  normal return is `stop()` on malformed `block` input).
- **C++ dependency:** **none.**

No other export, internal helper, or S3 method from `lmebayesCore` is
reachable from `lmbBlock()`.

---

## Look-alike dependencies (traced and ruled out)

These are easy to mistake for `lmebayesCore` dependencies because of naming,
`NAMESPACE` `import(lmebayesCore)`, or nearby documentation cross-references
— each was traced to its actual source and confirmed **not** to be
`lmebayesCore`:

| Symbol | Where it actually comes from | Why it looks like `lmebayesCore` |
|--------|-------------------------------|-----------------------------------|
| `.mrglmb_normalize_pfamily_lists()`, `.validate_pfamily_for_rlmb()` | **`glmbayesCore`** — `lmebayes/R/block_core_pfamily.R` binds them via `getFromNamespace(".mrglmb_normalize_pfamily_lists", "glmbayesCore")` (and same for the validator) | Same `multi_rlmb`/`mrglmb` naming family used elsewhere in the `lmebayesCore`/`glmbayesCore` mixed-model stack; also, `lmebayesCore` has its own **separate, independent** copy of `.mrglmb_normalize_pfamily_lists()` in `R/multi_rlmb.R` (used by *other* `lmebayesCore` entry points, e.g. `multi_rlmb()` there) — `lmbBlock()` does not call that copy. |
| `lmb()` | **`glmbayes`** (`glmbayes::lmb`) | `lmbBlock()`'s whole job is "one `lmb()` per block" — easy to assume `lmb()` itself routes through `lmebayesCore`/`glmbayesCore`. It does not: `glmbayes` has no `Imports`/`Depends`/`LinkingTo` on either package (confirmed in `glmbayes/DESCRIPTION`) and ships its own complete, independent compiled engine. |
| `lmebayes`'s `NAMESPACE` has `import(lmebayesCore)` (whole-namespace import) | N/A — this is a package-level import covering **all** of `lmebayes`, not evidence that `lmbBlock()` specifically uses more than one symbol from it. | Broad `import()` (vs. selective `importFrom()`) makes every `lmebayesCore` export unqualified-callable from anywhere in `lmebayes`, which can make grep-by-eye overestimate a given function's actual usage. |

---

## Related functions in the same file that are *not* part of `lmbBlock()`'s call graph

These live in `lmebayes/R/lmbBlock.R` alongside `lmbBlock()` and are
documented with `@seealso lmbBlock` / cross-linked from it, but are
**independent, user-invoked functions** — `lmbBlock()` does not call them:

| Function | Calls `lmebayesCore`? | Notes |
|----------|------------------------|-------|
| `block_check_identifiability_xy(x, block, ...)` | Yes — `lmebayesCore::normalize_block()` (same single function, `lmebayes/R/lmbBlock.R:403`) | Standalone preflight diagnostic for block/hyper identifiability; users call it separately, `lmbBlock()` never calls it. |
| `block_check_identifiability(formula, block, data, ...)` | No | Calls `.blmb_blocks_full_rank()` (internal), which uses only base R `qr()` — no `lmebayesCore`. |
| `.blmb_blocks_full_rank()`, `.blmb_blocks_full_rank_xy()` | No | Pure R (`stats::model.matrix`, `qr`). |

### Documentation-only cross-reference (never actually called)

`lmbBlock()`'s roxygen block contrasts it with
`\code{\link[lmebayesCore]{block_rNormalGLM}}` ("Gibbs conditional draws,
matrix API") purely as prose — explaining that `lmbBlock()` (independent
per-block `lm`-style fits) is a different modeling approach from
`lmebayesCore::block_rNormalGLM()` (a full conditional draw inside a Gibbs
sweep, matrix API, used by `lmerb()`/`glmerb()`'s Block~1 — see
`R_EXPORT_REACHABILITY.md`). `lmbBlock()` **never calls**
`block_rNormalGLM()` at runtime.

---

## Sibling function: `glmbBlock()`

Not asked for, but noted since it shares the same shared helper:
`glmbBlock()` (`lmebayes/R/glmbBlock.R`, the GLM-family counterpart to
`lmbBlock()`) calls the same `.blmb_formula_block_meta()` /
`.blmb_rows_to_data_subset()` helpers defined in `lmbBlock.R`, and therefore
has the **exact same** single `lmebayesCore` dependency:
`lmebayesCore::normalize_block()`. No additional `lmebayesCore` surface.

---

## Implication for the Stage 3 C++ dedup effort

Because `lmbBlock()` (and `glmbBlock()`) never reach any `lmebayesCore`
compiled code — only the pure-R `normalize_block()` — **neither function is
affected by the Stage 3 C++ bridge/dedup work** (see
`DESIGN_GLMBAYESCORE_BRIDGE.md`). Whatever happens to `lmebayesCore`'s C++
iid engine in Stages 3c-3f, `lmbBlock()`/`glmbBlock()` keep working
unchanged, since their only load-bearing dependency is a leaf R function with
no C++ underneath it.

---

## Verification method

- Read `lmebayes/R/lmbBlock.R` in full (`lmbBlock`, `print.blmb`, and every
  `.blmb_*` / `block_check_identifiability*` helper defined in that file).
- `grep -n "lmebayesCore" lmebayes/R/lmbBlock.R` → exactly 4 matches: 2 in
  roxygen `\code{\link[...]}` cross-references (prose only), 2 in actual code
  (`lmebayesCore::normalize_block(...)`, once inside `.blmb_formula_block_meta()`
  and once inside `block_check_identifiability_xy()` — the latter not on
  `lmbBlock()`'s call graph).
- Read `lmebayesCore/R/simfunction_block_utils.R` (`normalize_block()`'s
  full implementation) — confirmed pure R, no `.Call()`, no further
  `lmebayesCore`/`glmbayesCore` references.
- Read `lmebayes/R/block_core_pfamily.R` — confirmed
  `.mrglmb_normalize_pfamily_lists()` / `.validate_pfamily_for_rlmb()` are
  bound from `glmbayesCore` via `getFromNamespace()`, not `lmebayesCore`.
- Grepped `glmbayes/DESCRIPTION` and all of `glmbayes/R/`, `glmbayes/src/`
  for `lmebayesCore` → no matches; confirmed `glmbayes::lmb()` has no
  dependency on either core package.
- Grepped `lmebayes/NAMESPACE` for `lmebayesCore` → one `import(lmebayesCore)`
  (whole-namespace) plus 4 unrelated `importFrom()` entries
  (`build_mu_all`, `dGamma_list`, `glmerb_posterior_mode`,
  `lmerb_posterior_mean` — none used by `lmbBlock()`).

---

## Related files

| Topic | Path |
|-------|------|
| One-line summary this doc expands on | `inst/R_EXPORT_REACHABILITY.md` ("lmebayes — additional entry points" table) |
| `normalize_block()` implementation | `R/simfunction_block_utils.R` |
| `lmbBlock()` implementation (downstream package) | `lmebayes/R/lmbBlock.R` |
| `glmbBlock()` implementation (downstream package, sibling) | `lmebayes/R/glmbBlock.R` |
| pfamily-list validators actually from `glmbayesCore` | `lmebayes/R/block_core_pfamily.R` |
| Stage 3 C++ dedup bridge (why this doesn't affect `lmbBlock()`) | `inst/DESIGN_GLMBAYESCORE_BRIDGE.md` |

---

## Changelog

| Date | Note |
|------|------|
| 2026-07-20 | Initial assessment: `lmbBlock()`'s only `lmebayesCore` dependency is `normalize_block()` (pure R, no C++). Traced and ruled out `glmbayesCore` pfamily validators and `glmbayes::lmb()` as look-alike dependencies; documented `block_check_identifiability_xy()`/`glmbBlock()` as sharing the same single dependency without being on `lmbBlock()`'s own call graph. |
