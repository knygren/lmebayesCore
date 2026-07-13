# Per-group `dGamma()` measurement dispersion: the marginal Gamma and its bounds

How `dGamma_list()` (`R/dGamma_list_lmebayes_prior_setup.R`) builds a
per-group Independent Normal--Gamma (ING) prior on Block~1 measurement
dispersion `sigma^2_j`, and how it derives that prior's truncation window.
Companion to `inst/TAU2_ING_FORMULAS.md` (the same conjugate machinery for
`tau^2_k`) and to `lmebayes`'s `inst/GROUP_DISPERSION_HYPERPRIOR.md` /
`inst/ING_TRUNCATION_WINDOW.md`. Reproducible with
`data-raw/group_dgamma_bounds_derivation_check.R`.

---

## 0. Setup and notation

For each group (school) `j = 1..J`:

| Symbol | Code | Meaning |
|---|---|---|
| `n_j` | `g$n_j` | Observations in group `j` |
| `X_j` | `model.matrix(block_formula, dat_j)` | `n_j x p_re` design (RE-relevant predictors only; see `.lmebayes_block_formula_from_re()`) |
| `Y_j` | — | Response, group `j` |
| `bhat_j` | `coef(glm.fit(X_j, Y_j))` | Group `j`'s own OLS/WLS estimate |
| `mu_j` | `.lmebayes_block_formula_prior_mu()` | Prior mean: null-model intercept, zero slopes (defaults) |
| `V0_j` | `vcov(glm_full)` | Sampling covariance of `bhat_j` |
| `sd_tau` | `sqrt(diag(Sigma_ranef))` | Population RE sd, shared across groups |
| `w_j` | `pwt_measurement` (scalar or length-`J`) | Group `j`'s Block~1 prior weight |
| `n_prior,j` | `w_j/(1-w_j) * n_j` | Group `j`'s effective prior sample size |
| `n_combined,j` | `n_prior,j + n_j` | Total effective sample size |
| `p_re` | `length(design$re_coef_names)` | Number of RE components (2 in the worked example: intercept + `distracted_ppvt`) |

`.lmebayes_calibrate_ing_prior_measurement_group()` (`R/mixed_rmerb_helpers.R`)
builds, **independently for each group**, a coefficient-scale prior covariance

```
pwt_j    <- diag(V0_j) / (diag(V0_j) + sd_tau^2)          # per-coefficient shrinkage weight
Sigma_j  <- V0_j * outer(sqrt((1-pwt_j)/pwt_j), sqrt((1-pwt_j)/pwt_j))   # coefficient-scale prior cov
Sigma_0j <- Sigma_j / dispersion_classical_j               # dispersion-free (Zellner form)
```

and calls `compute_gaussian_prior(X_j, Y_j, weights=1, offset=0, dispersion=NULL,
n_effective=n_j, bhat=bhat_j, mu=mu_j, Sigma_0=Sigma_0j, Sigma=Sigma_j,
n_prior=n_prior,j, k=1)`.

---

## Part I -- The marginal Gamma distribution (per-group prior)

This section specializes `vignettes/Chapter-A12.Rmd` Sections 3.1, 3.3.1, and
3.3.4-3.3.5 (the general Normal--Gamma / ING / `dGamma()` calibration) to one
group's Block~1 design; see that vignette for the full derivation and proof
that this is the correct marginal (i.e. `beta` integrated out) Gamma law on
the residual precision.

**Marginal quadratic term** (Chapter A12 §3.1, specialized to group `j`,
weights `=1`):

```
S_marg,j = RSS_ols,j + (bhat_j - mu_j)' * Mi_j * (bhat_j - mu_j)
Mi_j     = (Sigma_0j + (X_j'X_j)^{-1})^{-1}
RSS_ols,j = sum((Y_j - X_j %*% bhat_j)^2)
```

The second term is a **penalty for how far the group's own OLS fit sits from
its prior mean `mu_j`** (relative to how confident the prior is, via
`Sigma_0j`), on top of the group's raw residual sum of squares.

**Calibrated point estimate** (`g$sigma2_hat`, `cal$dispersion`):

```
sigma2_hat,j = S_marg,j / (n_j - p_re)
```

Because `S_marg,j >= RSS_ols,j`, **`sigma2_hat,j >= sigma2_ols,j` whenever
`bhat_j != mu_j`** -- the quad-penalty term always inflates the calibrated
dispersion above the raw OLS residual variance. This gap is largest exactly
for groups whose own fit deviates most from the shared null-model prior mean
(see Part IV: this is also where `blup_infl,j` tends to be largest).

**The prior actually fed to the sampler** (`g$shape_ING`, `g$rate_gamma` --
this *is* `dGamma(shape = g$shape_ING, rate = g$rate_gamma, ...)`, group `j`'s
own Block~1 measurement-dispersion prior object):

```
shape_ING,j = (n_prior,j + 1)/2 + p_re/2                                    # k = 1
beta_star,j = (1 - w_j) * bhat_j + w_j * mu_j                               # Zellner blend at group j's own w_j
RSS_star,j  = sum((Y_j - X_j %*% beta_star,j)^2)
rate_gamma,j = 0.5 * RSS_star,j * (n_prior,j + p_re - 1) / (n_j - p_re)
```

Two things to note, both direct instances of Chapter A12 §3.3.5:

1. **This prior is calibrated per group**, using group `j`'s own `bhat_j`,
   `mu_j`, and `w_j` -- it is *not* one shared prior across all `J` groups
   (contrast with the "consistent design" `lmebayes`'s
   `inst/GROUP_DISPERSION_HYPERPRIOR.md` §8 argues for; see Part V).
2. `rate_gamma,j` uses `RSS_star,j` (RSS at the Zellner-blended coefficient),
   **not** `S_marg,j`. These differ whenever `w_j != 0`; `sigma2_hat,j`
   (used for the *bounds*, Part II) and `rate_gamma,j` (used for the *prior
   object*, fed to the sampler) are computed from genuinely different RSS
   quantities and must not be conflated.

`E[sigma^2_j] = rate_gamma,j / (shape_ING,j - 1)` under this prior alone
(before the sampler sees group `j`'s live data again) -- by design, `w_j`
(default `pwt_measurement = 0.01`) is small, so `n_prior,j << n_j` and this
prior carries little weight relative to what the Block~1 sampler will
subsequently learn from `(X_j, Y_j)`.

---

## Part II -- The truncation window (current `dGamma_list()` bounds)

### Why a window at all

Same three reasons as `lmebayes`'s `inst/ING_TRUNCATION_WINDOW.md` for
`tau^2_k`: a `disp_lower > 0` plug-in for convergence calibration, guaranteed
geometric ergodicity of the envelope accept/reject sampler, and a
sweep-independent (fixed) truncation so the invariant distribution does not
drift across Gibbs sweeps.

### Mean-matched construction

`dGamma_list.lmebayes_prior_setup()` (`R/dGamma_list_lmebayes_prior_setup.R`)
does **not** reuse `rate_gamma,j` for the window. Instead it builds a fresh
Gamma whose **mean is pinned to `sigma2_hat,j` by construction**, with
`n_combined,j` controlling only the spread:

```
shape_w,j <- (n_combined,j + 1)/2 + p_re/2
rate_w,j  <- sigma2_hat,j * (n_combined,j + p_re - 1)/2
```

Algebraically, `rate_w,j / (shape_w,j - 1) = sigma2_hat,j` **for every
`n_combined,j`** -- mean and spread are deliberately decoupled, so widening or
narrowing the window (by changing what feeds `n_combined,j`) never moves its
center.

### Asymmetric upper tail (`disp_upper_anchor = "blup"`, the default)

`lmer`'s BLUP for group `j`'s random effect can sit well away from group
`j`'s own OLS fit. `.lmebayes_group_blup_rss_inflation()` computes

```
blup_infl,j = RSS_blup,j / RSS_ols,j >= 1
```

(`RSS_blup,j` = residual sum of squares at the reference `lmer` fit's BLUP
coefficients for group `j`; `RSS_ols,j` at group `j`'s own OLS fit). The upper
tail is then anchored at an inflated rate:

```
rate_u,j <- rate_w,j * blup_infl,j        (rate_u,j == rate_w,j when disp_upper_anchor = "symmetric")

disp_lower,j <- 1 / qgamma(max_disp_perc,     shape = shape_w,j, rate = rate_w,j)
disp_upper,j <- 1 / qgamma(1 - max_disp_perc, shape = shape_w,j, rate = rate_u,j)
```

`.lmebayes_dgamma_window_cross_percentiles()` additionally reports cross
percentiles `R_lo`, `R_hi` diagnosing how asymmetric the resulting window is
(flagged via `asymmetric_R_lo` / `asymmetric_R_hi`, printed with
`warn_asymmetric = TRUE`).

### Worked numbers (39-school `big_word_club` fixture, `p_re = 2`, `max_disp_perc = 0.8`)

```
school 41: n_j=11  sigma2_hat=428.6  blup_infl=1.52  ->  disp_lower= 289.6  disp_upper= 823.3
school  9: n_j=11  sigma2_hat=137.1  blup_infl=1.26  ->  disp_lower=  92.6  disp_upper= 217.5
```

---

## Part III -- A proposed refinement: integrating over random-effect uncertainty

### Motivation

`lmer`'s BLUP shrinkage assumes **one pooled residual variance across all
groups**. When a group's *true* `sigma^2_j` differs substantially from the
pooled value, `lmer` mis-states that group's sampling uncertainty and
shrinks its BLUP too much (if `sigma^2_j` is smaller than pooled) or too
little (if larger) -- inflating or deflating `RSS_blup,j` relative to
`RSS_ols,j` for reasons that have nothing to do with what the *bounds* should
actually cover. `blup_infl,j` is a reasonable, cheap proxy for "how much
wider does the upper tail need to be", but it borrows its magnitude from a
model (`lmer`, homogeneous variance) that disagrees with the one being
fit (heterogeneous per-group variance) -- see Part IV for direct evidence
this heterogeneity is real.

### `EnvelopeCentering`'s trace correction (`src/EnvelopeCentering.cpp`)

At a fixed working dispersion `dispersion2`, the posterior mean of the random
effect `b_j` under prior `b_j ~ N(mu_j, Sigma_ranef)` and Gaussian likelihood
is the ridge/GLS estimator

```
b2_j = (X_j'X_j / dispersion2 + P)^{-1} (X_j'Y_j / dispersion2 + P %*% mu_j),   P = Sigma_ranef^{-1}
Cov(b2_j) = (X_j'X_j / dispersion2 + P)^{-1}
```

The **expected** residual sum of squares under this posterior -- not just the
RSS *at* the posterior mean -- follows from the standard bias--variance
(law-of-total-variance) decomposition of `E[(Y_j - X_j b_j)'(Y_j - X_j b_j)]`
over `b_j`'s posterior:

```
RSS_precomputed,j = ||Y_j - X_j %*% b2_j||^2  +  tr(X_j'X_j %*% Cov(b2_j))
                     ^^^^^^^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^^^^^^^^^^^^^^^^
                     RSS at the posterior mean    correction for b_j's own
                                                   remaining posterior uncertainty
```

`dispersion2` and `RSS_precomputed,j` are then updated to a fixed point (10
iterations in both `EnvelopeCentering.cpp` and the reproduction script):

```
dispersion2 <- RSS_precomputed,j / (n_j - p_re)     # data-only path -- see pitfall below
```

### Mean-matching the refined estimate

Exactly the Part II construction, with `sigma2_hat,j` replaced by the
converged `dispersion2,j`, and **no separate BLUP widening** -- the trace term
already prices in `b_j`'s uncertainty directly, symmetrically:

```
shape_w_prop,j <- (n_combined,j + 1)/2 + p_re/2      # unchanged
rate_w_prop,j  <- dispersion2,j * (n_combined,j + p_re - 1)/2   # mean-matched at dispersion2,j, not sigma2_hat,j

disp_lower_prop,j <- 1 / qgamma(max_disp_perc,     shape = shape_w_prop,j, rate = rate_w_prop,j)
disp_upper_prop,j <- 1 / qgamma(1 - max_disp_perc, shape = shape_w_prop,j, rate = rate_w_prop,j)
```

### The double-counting pitfall

**Do not** seed the `dispersion2` fixed point (or the final mean-match) with
group `j`'s own `shape_ING,j` / `rate_gamma,j` as if they were a fresh prior
to be updated with `n_j/2`, `RSS_precomputed,j/2` on top. `shape_ING,j` /
`rate_gamma,j` (Part I) are **already** calibrated from group `j`'s full `n_j`
observations (via `RSS_star,j`, itself derived from `bhat_j`) -- just scaled
down to an equivalent `n_prior,j`-sized weight for a different purpose (the
prior object handed to the sampler). Running the same `n_j` observations
through the update a second time inflates the effective weight of the
group's own data well beyond what `n_prior,j` is meant to represent, and
produces artificially tight windows. Concretely, on the 39-school fixture this
mistake pushed school 9's BLUP point estimate (`sigma2_blup = 160.4`) outside
its own proposed upper bound (`disp_upper_prop = 157.2`) -- a symptom that a
group's own reference-model point estimate should essentially never land
outside an 80%-mass window built for exactly that group. Mean-matching
directly at `dispersion2,j` with the *same* `n_combined,j` convention Part II
already uses (rather than re-deriving shape/rate through `shape_ING,j` +
`n_j/2`) avoids this: `n_j` enters through `RSS_precomputed,j`'s point
estimate exactly once.

### Worked numbers (same fixture and setting)

| group | n_j | sigma2_ols | sigma2_blup | blup_infl | dispersion2 (proposed center) | disp_lower (current) | disp_upper (current) | disp_lower (proposed) | disp_upper (proposed) | %delta lower | %delta upper |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 41 | 11 | 283.6 | 432.3 | 1.52 | 472.6 | 289.6 | 823.3 | 319.3 | 595.4 | +10.3% | -27.7% |
| 33 | 12 | 209.4 | 292.5 | 1.40 | 339.7 | 225.3 | 573.7 | 233.5 | 425.8 | +3.7% | -25.8% |
| 6  | 11 | 272.9 | 373.3 | 1.37 | 404.4 | 263.5 | 672.2 | 273.3 | 509.6 | +3.7% | -24.2% |
| 9  | 11 | 127.4 | 160.4 | 1.26 | 156.6 | 92.6 | 217.5 | 105.8 | 197.3 | +14.2% | -9.3% |

With the pitfall fixed, the BLUP point estimate falls inside the proposed
window for **all 39 groups** (verified by
`data-raw/group_dgamma_bounds_derivation_check.R`). Across all 39 groups,
median `%delta lower = +11.2%` (proposed lower bound is looser than current),
median `%delta upper = +6.0%`, but the three highest-`blup_infl` groups (6,
33, 41) get **narrower** upper bounds than the current BLUP-inflated ones
(-24% to -28%) while most other groups get modestly wider ones (up to +23%).

---

## Part IV -- Empirical validation against an independent classical estimator

Both `sigma2_hat,j` (Part I's quad-penalized point estimate) and the
`blup_infl,j`-driven asymmetry (Part II) rest on the premise that residual
variance genuinely differs across groups -- `lmer` assumes otherwise. This is
directly testable: fit the *same* random-slope structure with `nlme::lme()`,
adding a per-group residual-variance term via `varIdent()`:

```r
fit_hetero <- nlme::lme(
  score_ppvt ~ 1 + distracted_ppvt,
  random  = list(school_id = nlme::pdDiag(~ 1 + distracted_ppvt)),
  weights = nlme::varIdent(form = ~ 1 | school_id),
  data    = dat
)
```

**Heteroscedasticity is statistically supported.** LRT of `lme4` (pooled
`sigma^2`) vs. `nlme` (per-group `sigma^2`), same RE structure: `57.36` on
`38` df, `p = 0.023`.

**`nlme`'s REML per-group variances track both our point estimates closely**
(`cor(sigma2_ols, sigma2_nlme) = 0.985`, `cor(sigma2_hat, sigma2_nlme) =
0.995`), confirming the per-group heterogeneity direction and rough magnitude
independently of this package's own calibration.

**The quad penalty (`sigma2_hat,j > sigma2_ols,j`) is validated exactly where
it is largest.** `nlme`'s estimate falls strictly between `sigma2_ols,j` and
`sigma2_hat,j` for 16 of 39 groups; restricting to the 11 groups with the
largest `blup_infl,j` (equivalently, the largest quad-penalty gap
`sigma2_hat,j - sigma2_ols,j`), it holds for **all 11** (correlation between
`blup_infl,j` and "nlme falls between" = `0.69`; between the quad-penalty gap
and "nlme falls between" = `0.73`). For groups with `blup_infl,j ~ 1` (quad
penalty negligible), `nlme` instead lands slightly below both -- its own
REML shrinkage-toward-pooled effect on noisy small-`n_j` variance estimates,
unrelated to the quad penalty.

---

## Part V -- Status and open questions

- **Part III is not implemented in `dGamma_list()`.** This document records
  the derivation and validation for a future change; `disp_upper_anchor =
  "blup"` (Part II) remains the shipped default.
- **Tension with `lmebayes`'s `inst/GROUP_DISPERSION_HYPERPRIOR.md` §8.**
  That note argues the "consistent design" keeps **one shared prior** across
  all groups, varying only the truncation window per group, to avoid
  double-counting each group's data in its own prior. `dGamma_list()` instead
  calibrates `shape_ING,j` / `rate_gamma,j` **per group** (Part I). The
  default `pwt_measurement = 0.01` keeps `n_prior,j << n_j`, so the resulting
  double-count in the *prior itself* is proportionally small -- but it is not
  zero, and is a distinct issue from the Part III bounds double-count (which
  is a full, not partial, double-count and must be avoided entirely).
- Part III's `n_combined,j`-based shape convention was carried over unchanged
  from Part II for direct comparability; whether it remains the right spread
  choice once the center is `dispersion2,j` rather than `sigma2_hat,j` is an
  open question.

---

## References

- `vignettes/Chapter-A12.Rmd` -- the general Normal--Gamma / ING / `dGamma()`
  marginal-Gamma derivation this document specializes (Part I).
- `inst/TAU2_ING_FORMULAS.md` -- the same conjugate Gamma--Normal update
  mechanics applied to `tau^2_k` instead of measurement dispersion.
- `lmebayes`'s `inst/ING_TRUNCATION_WINDOW.md` -- the truncation-window
  rationale this document's Part II specializes to per-group bounds.
- `lmebayes`'s `inst/GROUP_DISPERSION_HYPERPRIOR.md` -- precursor exploration
  of group-level dispersion heterogeneity and the double-counting hazard
  (§8), on the same `big_word_club` fixture.
- Pinheiro, J.C. and Bates, D.M. (2000). *Mixed-Effects Models in S and
  S-PLUS.* Springer. -- `varIdent()` heteroscedastic variance structures used
  in Part IV.

## Runnable check

```r
Rscript data-raw/group_dgamma_bounds_derivation_check.R
```

Reproduces the Part II/III worked-numbers table (all 39 groups) and the Part
IV `nlme` validation statistics from the live `big_word_club` fixture.
