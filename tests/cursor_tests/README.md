# cursor_tests — manual regression via **lmebayes** demos

This folder is **not** part of `R CMD check`. Do **not** add `test-*.R`, smoke
scripts, or ad-hoc fixtures here.

## Policy

All manual regression for the six matrix-level **reg route** engines in
**glmbayesCore** must use the **lmebayes** demos and help examples listed in
[lmebayes/README.md — Examples and Demos](https://github.com/knygren/lmebayes/blob/main/README.md)
(plus `Ex_23` below for the full LMM 2×2). Workflow everywhere:

`model_setup()` → `Prior_Setup_lmebayes()` → `pfamily_list(ps)` →
`lmerb()` / `glmerb()` (which call `rlmerb()` / `rglmerb()` →
`.lmebayes_run_*_engine()` → `REG_ROUTE_TABLE`).

Requires **bayesrules** (and usually **lme4**) where noted in each demo.

## Six route engines ↔ demos

Formula drivers dispatch by `(family, σ² mode, Block~2 ING)`. Matrix exports
are reached only through **`rlmerb()`** / **`rglmerb()`** unless you call them
directly with the same **`prior_list`** / **`pfamily_list`** objects built from
the demo’s **`Prior_Setup_lmebayes()`** output.

| Matrix route (`glmbayesCore`) | Route key | σ² | Block~2 | Formula driver | Run this demo / example |
|-------------------------------|-----------|----|---------|----------------|-------------------------|
| `rLMMNormal_reg_known_vcov()` | `lmm_fixed_known` | fixed scalar | all `dNormal` | `lmerb()` / `rlmerb()` | `example("lmerb")`; `demo("Ex_14_lmerb_Sleepstudy", package = "lmebayes")`; `demo("Ex_12_lmerb_BigWordClub", package = "lmebayes")`; **Ex_23 case 1** |
| `rLMMNormal_reg_estimated_vcov()` | `lmm_fixed_estimated` | fixed scalar | ≥1 ING | `lmerb()` / `rlmerb()` | `demo("Ex_20_lmerb_ING_pilot", package = "lmebayes")`; `demo("Ex_21_lmerb_ING_BigWordClub", package = "lmebayes")`; **Ex_23 case 2** |
| `rLMMindepNormalGamma_reg_known_vcov()` | `lmm_gamma_known` | `dGamma()` Block~1 | all `dNormal` | `lmerb()` / `rlmerb()` | **`demo("Ex_23_lmerb_joint_posterior_mode_four_cases", package = "lmebayes")` case 3** (ICM); for Gibbs, re-run case 3 with `simulate = TRUE`, small `n` |
| `rLMMindepNormalGamma_reg_estimated_vcov()` | `lmm_gamma_estimated` | `dGamma()` Block~1 | ING | `lmerb()` / `rlmerb()` | **Ex_23 case 4** (ICM); for Gibbs, re-run case 4 with `simulate = TRUE`, small `n` |
| `rGLMM_reg_known_vcov()` | `glmm_known` | none (Poisson/binomial) | all `dNormal` | `glmerb()` / `rglmerb()` | `example("glmerb")`; `demo("Ex_14_glmerb_airbnb_small", package = "lmebayes")`; `demo("Ex_16_glmerb_book_banning", package = "lmebayes")` |
| `rGLMM_reg_estimated_vcov()` | `glmm_estimated` | none | ING | `glmerb()` / `rglmerb()` | `demo("Ex_22_glmerb_book_banning_ING", package = "lmebayes")` |

**Ex_23** = `demo("Ex_23_lmerb_joint_posterior_mode_four_cases", package = "lmebayes")`
— documents all four Gaussian LMM subcases (fixed vs dGamma σ² × dNormal vs ING
Block~2). It is the only shipped demo for the two **dGamma σ²** LMM routes; it
uses `simulate = FALSE` (ICM / joint mode). Use the same calls with
`simulate = TRUE` when testing the ING Block~1 sweep engines end-to-end.

## Full demo index (from **lmebayes** README)

Help examples (ICM / setup; fast):

```r
example("lmerb")    # big_word_club Gaussian LMM
example("glmerb")   # airbnb_small Poisson GLMM
```

Gibbs demos (stored draws; may take minutes):

```r
demo("Ex_14_lmerb_Sleepstudy", package = "lmebayes")
demo("Ex_14_glmerb_airbnb_small", package = "lmebayes")
demo("Ex_12_lmerb_BigWordClub", package = "lmebayes")
demo("Ex_13_glmerb_Airbnb", package = "lmebayes")
demo("Ex_16_glmerb_book_banning", package = "lmebayes")
demo("Ex_19_glmerb_book_banning_state_covariates", package = "lmebayes")
demo("Ex_20_lmerb_ING_pilot", package = "lmebayes")
demo("Ex_21_lmerb_ING_BigWordClub", package = "lmebayes")
demo("Ex_22_glmerb_book_banning_ING", package = "lmebayes")
demo("Ex_23_lmerb_joint_posterior_mode_four_cases", package = "lmebayes")
```

## Related **glmbayesCore** checks

- **`tests/testthat/`** — package check suite only (`test-prior-setup-poisson-conj.R` today).
- **`inst/examples/`**, **`demo/`** — iid GLM / block utilities; not reg-route coverage.
- **`data-raw/`** — maintainer scratch; not regression tests.
