# Role of examples, vignettes, manual tests, demos, and tests

This note describes how **lmebayesCore** splits work across five layers.
Each layer has a different audience, runtime budget, and CRAN obligation.
They are complementary: none replaces another.

| Layer | Runs on CRAN check? | Typical `n` / runtime | Primary purpose |
|-------|---------------------|------------------------|-----------------|
| **Examples** (`man/*.Rd`, `@examples`) | Yes — **each** example must finish within **5 seconds** | Small (`n` in the tens–low hundreds; ICM-only where possible) | Show that exported API calls run and return sensible structure |
| **Vignettes** (`vignettes/*.Rmd`) | Yes — but **heavy computation should be pre-generated** | Large when illustrating theory (`n = 10000` in stored artifacts) | Teach theory, workflows, and interpretation with stable output |
| **Manual tests** (`data-raw/`, `tests/cursor_tests/`, maintainer scripts) | **No** | As long as needed; run selectively | Verify theory, numerical accuracy, and full-route behaviour |
| **Demos** (`demo/`) | No (user runs `demo()`) | Medium (`n ≈ 1000` or less) | Let users explore models interactively without vignette-length waits |
| **Tests** (`tests/testthat/`) | Yes | Short (`n` small; tight tolerances on *rates*, not full posteriors) | Guard regressions on invariants: acceptance, convergence calibration, API contracts |

---

## 1. Examples — short, CRAN-bound

**Policy:** every `\examples` block must be safe to run during `R CMD check`.
CRAN enforces a **5-second ceiling per example** (each `\examples` block
run separately during check). Treat that as a hard budget for every help
page, not an average across the package.

**Do:**

- Prefer **ICM / joint-mode** paths (`simulate = FALSE`, or engines that
  return a mode without storing `n` draws) when the example’s goal is setup
  or a single posterior summary.
- Keep **`n` small** when draws are unavoidable (enough to show dimensions and
  class, not to estimate tail probabilities precisely).
- Use `\donttest{}` only when the example is genuinely slow but still
  runnable; do not use it to hide broken code. Prefer making the example
  fast.
- Gate progress noise with existing options (`progbar = FALSE`,
  `verbose = FALSE`) so check output stays clean (see
  `.cursor/rules/cran-rd-value-and-console.mdc`).

**Do not:**

- Run multi-hour Gibbs jobs, `n = 10000` sweep replication, or large
  benchmark grids inside `@examples`.
- Duplicate vignette-length narrative in help pages.

**Role:** examples answer “does this function run and return the documented
structure?” — not “does the sampler match Nygren (2020) at `n = 10000`?”

---

## 2. Vignettes — long when pre-generated

**Policy:** vignettes may be **as long and detailed as the subject requires**,
provided **`R CMD check` does not re-run the expensive parts**.

**Pattern used in this package:**

1. Maintainer script under `data-raw/make_Chapter-*.R` (or similar) runs the
   heavy fit once — e.g. `n = 10000`, full sweep history, reference
   `two_block_rate` objects.
2. Output is stored in `inst/extdata/*.rds` (summaries, sweep history,
   convergence metadata, stripped design objects for plots).
3. The vignette **loads and displays** the artifact; chunks that only knit
   tables, theory, and plots stay fast.

Example: `vignettes/Chapter-C02.Rmd` and `Chapter-C04.Rmd` share the same
Gaussian design; C02 uses the exact iid route (~seconds at `n = 10000`), while
C04 stores a two-block Gibbs run that took hours to produce. Both vignettes
render in check time because the sampler is not re-executed during knit.

**Do:**

- Reproduce **theoretical properties** (TV bounds, variance ratios, tail
  probabilities, sweep budgets) with **`n` large enough** that Monte Carlo
  noise does not dominate the story.
- Document the generation script and seed in the vignette (“produced by
  `data-raw/make_Chapter-C04.R` with `n = 10000`, `set.seed(1)`”).
- Keep vignettes **model- or theory-focused** (C-series for Gaussian theory;
  B-series for engine documentation).

**Do not:**

- Leave `\eval=TRUE` chunks that call `rlmerb(n = 10000, …)` with
  `TWO_BLOCK_GIBBS` unless you have measured knit time and it fits check
  policy on all CRAN platforms.
- Patch generated `man/*.Rd` or HTML previews without updating the `.Rmd`
  source.

**Role:** vignettes answer “what does this package mean, and what should I
expect from a correct run?” — with publication-quality numbers.

---

## 3. Manual tests — long, selective, off-CRAN

**Policy:** manual tests are **maintainer-owned** regression and theory
checks. They are **not** part of `R CMD check` and may run for hours.

**Locations:**

| Location | Use |
|----------|-----|
| `data-raw/make_*.R`, `data-raw/_scratch_*.R` | Precompute vignette artifacts; one-off probes; timing studies |
| `tests/cursor_tests/` | Documented manual regression paths (see `README.md` there); **not** run by `test_check()` |
| Ad hoc scripts at repo root / `data-raw/` | Diagnostics during development — delete or promote deliberately |

**Do:**

- Re-run manual tests **selectively** when changing ergodicity code, sweep
  calibration, `plot_var_convergence()`, or artifact generation.
- Compare **full numerical output** against theory: eigenvalues, certified
  `m_convergence`, TV bound tables, exact vs empirical covariance ratios.
- Use the **same models** as vignettes but allow **longer** runs or extra
  variants (different seeds, tolerances, pilot paths).

**Do not:**

- Add new files under `tests/testthat/` without explicit approval (see
  `.cursor/rules/testthat-and-scratch-code.mdc`).
- Treat `data-raw/_scratch_*.R` as shipped regression tests — they are
  disposable unless promoted to `make_*.R` or documented manual tests.

**Role:** manual tests answer “are we still correct, including edge cases
CRAN will never exercise?” — the ground truth for vignette artifacts.

---

## 4. Demos — shortened vignettes for users

**Policy:** demos are **interactive, user-invoked** scripts (`demo("…")`).
They should let someone **try a model end-to-end** without committing to a
vignette-length or manual-test-length run.

**Guidelines:**

- Target **`n ≈ 1000`** (or lower) for Gibbs / multi-sweep routes unless the
  demo is ICM-only.
- Mirror the **same workflow** as the corresponding vignette or manual test:
  `model_setup()` → `Prior_Setup_GLMM()` → `pfamily_list()` →
  `lmerb()` / `glmerb()` / matrix route — but with fewer draws and less
  printed detail.
- Keep **`demo/00Index`** and header comments aligned with the five LMM/GLMM
  route matrix (see `tests/cursor_tests/README.md` for the reg-route map).

**Relationship to vignettes:**

| Vignette / manual test | Demo |
|------------------------|------|
| `n = 10000`, stored artifact, full theory | `n = 1000`, live run, core output only |
| Explains *why* (bounds, plots, appendices) | Shows *how* (callable script) |

**Role:** demos answer “can I run this model myself in a few minutes?” —
exploration, not certification.

---

## 5. Tests (`testthat`) — short invariants with optional anchors to manual work

**Policy:** `tests/testthat/` runs on **every** check. Tests must be
**fast**, **deterministic**, and focused on **properties that must not
regress** — not on reproducing full vignette numerics.

**Good test targets:**

- **API contracts:** argument validation, return classes, NAMESPACE-visible
  behaviour, error messages.
- **Convergence machinery:** `two_block_rate()` eigenvalues match a small
  fixed fixture; `two_block_tv_bound()` monotone in `l`; certified
  `m_convergence` formula (+1 half-step) on a toy precision matrix.
- **Acceptance / mixing proxies:** where MCMC is involved, use **tiny `n`**
  and check **rates or bounds** (e.g. sweep count within a band, variance
  ratio trending toward 1 by sweep `m`, envelope acceptance in `[0,1]`) —
  not that posterior means match a stored `n = 10000` table to three decimals.

**Anchoring to manual tests (optional but useful):**

When a manual run establishes a reference (e.g. `lambda_star = 0.9308`,
Theorem 3 certifies 24 sweeps at `tv_tol = 0.01`), a **unit test** can check
a **cheap surrogate** on the same design:

- same **`lambda_star`** on a precomputed rate object or small fixture;
- **`m_convergence`** within ±1 of the manual value;
- variance-ratio at final sweep **below 1** and **above 0.9** at `n = 50`
  chains — not identical to the manual `n = 10000` plot.

Avoid brittle tests that fail on unrelated platform noise; prefer inequalities
and structural checks derived from theory.

**Do not:**

- Run `n = 1000` Gibbs replication suites in `testthat` by default.
- Duplicate entire vignette knit paths inside tests.

**Role:** tests answer “did we break a invariant the package relies on?” —
continuous guardrails, not a second manual test suite.

---

## Workflow summary

```
Theory / full accuracy          Manual test (data-raw, cursor_tests)
        │                                    │
        │  precompute                        │  selective re-run
        ▼                                    ▼
   inst/extdata/*.rds  ──────────►  Vignette (load artifact, teach)
        │
        │  shorten n, simplify output
        ▼
      Demo (user exploration)
        │
        │  extract minimal invariant
        ▼
   testthat (fast regression)

Examples: parallel minimal API smoke — each must stay within CRAN's 5 s limit
```

**When changing sampler or ergodicity code:**

1. Update **manual test** / regeneration script; confirm theory numbers.
2. Refresh **`inst/extdata`** if vignette output shifts.
3. Adjust **vignette text** if interpretation changes.
4. Update **demo** only if user-facing workflow changed.
5. Tighten or add **testthat** checks for new invariants — not full posteriors.
6. Keep **examples** fast; extend `\examples` only if the new API needs a
   minimal runnable call.

---

## Related policy in this repo

- `.cursor/rules/testthat-and-scratch-code.mdc` — where scratch and test
  files may live
- `.cursor/rules/cran-rd-value-and-console.mdc` — `\value`, `message()`
  vs `print()` on CRAN
- `tests/cursor_tests/README.md` — manual reg-route demo matrix
- `data-raw/make_Chapter-C*.R` — C-series artifact generation (vignette
  inputs)
