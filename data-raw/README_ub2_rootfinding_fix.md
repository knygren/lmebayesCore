# Port to `glmbayes`: exact (root-finding) `UB2_Min_j` instead of endpoint-only

**Status: IMPLEMENTED and validated in glmbayesCore, AND ported to
`glmbayes`.** This document records what changed so the same patch can be
copied to `glmbayes`'s (near-identical) copy of
`src/EnvelopeDispersionBuild.cpp`. It **replaces and supersedes**
`data-raw/README_glmbayes_g_prior_guard.md` (that guard has been removed from
glmbayesCore now that the actual root cause is fixed -- see "Do NOT port the
Zellner g-prior guard" below).

## Follow-up: near-isotropic fast path

**Status: IMPLEMENTED and validated in glmbayesCore; ported to `glmbayes`
(code change applied -- not yet compiled/tested there; see "Porting status"
below).**

Per Claim 7 part 3 (isotropic case) combined with parts 1/2 (any critical
point `t* < lambda_max(K)`; any inflection point `t* >= 2*lambda_min(K)`),
whenever `kappa(K) = lambda_max(K)/lambda_min(K) <= 2`, these two ranges
cannot overlap, so no interior local minimum of `UB2_j(d)` is possible and
the plain endpoint comparison is *exactly* the true minimum -- not just a
heuristic. `bound_ub2_over_dispersion()` now computes `kappa(K)` once (right
after the shared `K`/`K_eigval` eigendecomposition, before the per-face
loop -- `kappa(K)` does not depend on the face `j`) and, when
`kappa(K) <= 2` (with a `1e-8` relative cushion for eigendecomposition
noise), skips `ub2_exact_detail::ub2_min_exact_1d()` entirely for every
face, falling straight through to the (now provably exact) endpoint-only
comparison. This avoids the per-face root-finding cost whenever the
coefficient prior is exactly or nearly a Zellner g-prior for the
design/weights in use, which is expected to be the common case.

Both `kappa_K` and `K_is_near_isotropic` are surfaced in the `diagnostics`
list returned by `EnvelopeDispersionBuild()` (and therefore in
`rindepNormalGamma_reg_with_envelope()$diagnostics`) purely for
inspection/testing; they are not consumed by the sampler itself. See
`data-raw/validate_near_isotropic_fastpath.R` for a validation script
confirming: (1) an exact Zellner g-prior gives `kappa_K == 1` and triggers
the fast path; (2) a mild, non-Zellner vector-`pwt`-style perturbation with
`1 < kappa_K <= 2` still triggers the fast path; (3) a strongly anisotropic
prior (`kappa_K` in the hundreds of thousands) does not, and still runs the
full exact search.

**Ported to `glmbayes` (`C:\Rpackages\glmbayes\src\EnvelopeDispersionBuild.cpp`):**
added the same `kappa_K`/`K_is_near_isotropic` computation right after
`glmbayes`'s copy of the `K`/`K_eigval` eigendecomposition in its
`bound_ub2_over_dispersion()`, and gated its per-face
`ub2_exact_detail::ub2_min_exact_1d()` call on `!K_is_near_isotropic` the same
way. Also surfaced both diagnostics through `bound_ub2_over_dispersion()`'s
return list and into `EnvelopeDispersionBuild()`'s `diagnostics` list (the
glmbayesCore-only convenience noted above -- ported anyway for parity/
testability), plus a matching verbose `Rcout` line reporting `kappa(K)` and
`K_is_near_isotropic` after the `bound_ub2_over_dispersion()` call. Also
added a `NEWS.md` entry under "Independent Normal-Gamma simulation"
("Updates to independent normal gamma simulation to better handle
non-isotropic priors with highly differentiated implied pweights across
dimensions"). **Not yet compiled or tested in `glmbayes`** -- the user will
build and run checks separately; re-run (or adapt)
`data-raw/validate_near_isotropic_fastpath.R` against `glmbayes` before
relying on this.

### Practical design guidance: orthogonal design + close `pwt`'s is a free fast path

A closed-form corollary of `K_is_near_isotropic`, worked out analytically
(not by search) for the case `Q = X^TWX` diagonal (an orthogonal design --
either by actual experimental design, or after orthogonalizing/decorrelating
predictors as a preprocessing step -- with *arbitrary*, possibly unequal,
per-coordinate scales/norms), combined with a diagonal coefficient prior
(`Prior_Setup`'s vector-`pwt` mechanism, `Sigma = V0 \odot outer(s,s)`,
`s_i = sqrt((1-pwt_i)/pwt_i)`):

- `K = Q^{-1/2} P_coef Q^{-1/2}` is then exactly diagonal, with
  `lambda_i(K) = (1/sigma_hat^2) * pwt_i/(1-pwt_i)` -- the per-coordinate
  design scale `q_i` cancels out exactly (because `Prior_Setup`'s vector-`pwt`
  construction is *defined* relative to each coordinate's own
  likelihood-implied scale, `V0_ii \propto 1/q_i`, which is precisely what
  makes `pwt_i` a scale-free "prior weight" in the first place).
- Because `x/(1-x)` is monotonic, this ordering-preserving map means
  `lambda_max(K)` and `lambda_min(K)` come *only* from the single largest and
  single smallest `pwt_i` in the whole vector -- every intermediate `pwt_i`
  is irrelevant, for any number of coefficients `p`.
- So `kappa(K) <= 2` (the fast-path certificate) reduces to a simple,
  checkable, dimension-free, scale-free condition on just the two extreme
  `pwt`'s:
  ```
  pwt_max/(1-pwt_max) <= 2 * pwt_min/(1-pwt_min)
  ```
  which for small `pwt`'s (the common case: weak-ish priors relative to the
  data) simplifies further to the memorable rule of thumb
  `pwt_max <~ 2 * pwt_min` (slightly conservative; the exact bound is a
  touch tighter, by a factor `1/(1+pwt_min)`).

**Practical takeaway:** if a user (a) arranges for the design to be
orthogonal -- via actual experimental design, or by decorrelating
predictors -- and (b) keeps a diagonal/vector-`pwt` prior's most-informed and
least-informed coefficients within roughly a factor of 2 of each other in
`pwt`-odds, the near-isotropic fast path is *guaranteed*, regardless of `p`
and regardless of the design's raw per-coordinate scales. Caveats:

- This is a *sufficient*, not necessary, condition for safety in the
  correctness sense -- violating it only costs the (now-correct, exact)
  root-finding search, not correctness; the earlier `UB2_Min` fix already
  makes every case correct regardless of `kappa(K)`. This is purely a
  performance/design guideline.
- "Orthogonal" means orthogonal in the *weighted* inner product `X^TWX` at
  the weights actually used by the sampler; for the fixed-weight Gaussian/ING
  case this holds throughout sampling, but would need re-checking for
  IRLS-type weights that change with the fit (other GLM families).
- Orthogonalizing/decorrelating predictors changes what each coefficient
  (and hence each `pwt_i`) means -- `pwt_i` then applies to a rotated
  combination of the original variables, not the named coefficient a user
  may hold an actual prior belief about. This is a real trade-off for
  observational data; a genuinely designed experiment with orthogonal
  contrasts chosen up front doesn't have this problem, since orthogonality
  and interpretability coincide by construction.
- The moment `Q` has *any* off-diagonal correlation, this clean
  scale-cancellation breaks down and `kappa(K)` can grow far faster than the
  raw `pwt` spread would suggest. Writing `t_i = sqrt(pwt_i/(1-pwt_i))` and
  `tau = t_max/t_min` (so `tau^2` is the `pwt`-odds ratio -- `tau^2 <= 2` is
  exactly the diagonal-`Q` near-isotropic rule above), a worst-case
  operator-norm argument gives the general (correlated-`Q`) sufficient
  condition
  ```
  (tau - 1) * sqrt(kappa(Q)) <= 3 - 2*sqrt(2) ~= 0.1716  =>  kappa(K) <= 2
  ```
  i.e. the *tolerance on how different the `pwt`'s can be* shrinks as
  `1/sqrt(kappa(Q))` -- the design's own collinearity. Concretely: `tau = 2`
  (a `pwt`-odds ratio of `4`) is harmless for a well-conditioned design
  (`kappa(Q) = 3`, e.g. correlation `0.5` between two predictors) but was
  shown earlier in this design discussion to push `kappa(K)` to about `19`
  once the predictors are strongly correlated (`kappa(Q) = 19`, e.g.
  correlation `0.9`) -- well past the `<=2` safety bar. So the clean
  `pwt`-only rule above is strictly an *orthogonal-design* result; the moment
  predictors are correlated, `pwt`-closeness must be judged jointly with
  `kappa(Q)`, not on its own.

### Unifying name: "near-Zellner g-prior"

The classical Zellner g-prior is `Sigma = g*sigma_hat^2*Q^{-1}` -- a single
scalar `g` (equivalently, a single shared `pwt`) multiplying the inverse
Fisher information, i.e. `T = t*I` exactly. That gives `K = t^2*I` exactly,
`kappa(K) = 1`, for *any* design `Q`, not just orthogonal ones -- it is the
special case the pre-erratum Claim 7 part 3 implicitly assumed.

The `kappa(K) <= 2` near-isotropic certificate is exactly the natural
generalization of this: the set of priors whose `T` is a **bounded relative
perturbation away from scalar**, `T = t_min*(I + E)`, with the perturbation
budget `‖E‖` shrinking as `1/sqrt(kappa(Q))`. It is reasonable to call this
whole class **"near-Zellner g-priors"** -- priors that are a `T`-perturbation
of a true Zellner g-prior, where "how close is close enough" is design-aware
rather than a fixed tolerance on `pwt` alone:

- For an orthogonal design (`kappa(Q) = 1`), the near-Zellner neighborhood is
  at its widest: any diagonal/vector-`pwt` prior with
  `pwt_max/(1-pwt_max) <= 2*pwt_min/(1-pwt_min)` qualifies (the "Practical
  design guidance" rule above).
- As `kappa(Q)` grows, the neighborhood shrinks like `1/sqrt(kappa(Q))` and,
  in the limit of an ill-conditioned design, collapses to *only* the exact
  scalar-`pwt` Zellner prior itself (`tau -> 1`).
- Critically, "near-Zellner in `pwt`-space" (small `tau`) does **not** imply
  "near-Zellner in the sense that matters" (`kappa(K) <= 2`) unless the
  design is also reasonably well-conditioned -- the `tau=2`/`kappa(Q)=19`
  counterexample above (`kappa(K) ≈ 19`) is the cautionary case: a `pwt`
  spread that looks like a trivial perturbation from a shared-`pwt` Zellner
  prior can still land far outside the fast-path-safe region once the design
  is collinear enough.

## Porting status (glmbayes, `C:\Rpackages\glmbayes`)

Done:
- `src/EnvelopeDispersionBuild.cpp`: includes + `ub2_exact_detail` namespace
  + rewritten `bound_ub2_over_dispersion` ported verbatim (adapted only to
  keep glmbayes's existing `disp_min_ub2`/`ub2_min`-only return shape, i.e.
  no `ub2_at_low`/`ub2_at_upp` fields, since glmbayes never had those
  diagnostics). Compiled successfully via `R CMD INSTALL`.
- `vignettes/Chapter-A07.Rmd`: Remark 5.5.7 and Claim 7 erratum + corrected
  proofs, plus the new `5.5.5'`/Remark 5.5.9 section, ported verbatim
  (glmbayes ships its own copy of this vignette).
- `tests/testthat/test-lmb-non-zellner.R`: extended with a new
  `test_that("lmb: strongly anisotropic Independent Normal-Gamma prior does
  not trigger UB2 sign violations", ...)` block -- 10 reps of
  `lmb()` with a `diag(2000, 0.05)` coefficient prior (condition number
  40,000) and `n = 14`, asserting no error is thrown and all posterior
  means are finite. Full `testthat::test_dir()` run: all tests pass (OpenCL
  tests skip "On CRAN" as expected).
- No Zellner g-prior guard existed in `glmbayes` to begin with (confirmed by
  grep) -- nothing to remove there.

Done (near-isotropic fast path, applied after the above; **not yet
compiled/tested** -- left for the user to build and check):
- `src/EnvelopeDispersionBuild.cpp`: `kappa_K`/`K_is_near_isotropic`
  computation added right after `bound_ub2_over_dispersion`'s `K`/`K_eigval`
  eigendecomposition; per-face `ub2_exact_detail::ub2_min_exact_1d()` call
  gated on `!K_is_near_isotropic`; both diagnostics added to
  `bound_ub2_over_dispersion()`'s and `EnvelopeDispersionBuild()`'s returned
  lists; matching verbose `Rcout` line added.
- `NEWS.md`: added a short entry under a new "Independent Normal-Gamma
  simulation" heading in the `0.9.6.9000 (development)` section: "Updates to
  independent normal gamma simulation to better handle non-isotropic priors
  with highly differentiated implied pweights across dimensions."
- Not yet done: compiling/installing `glmbayes` with this change, and
  re-running/adapting `data-raw/validate_near_isotropic_fastpath.R` and the
  existing anisotropic-prior test in `tests/testthat/test-lmb-non-zellner.R`
  against it.

## Problem (recap)

Chapter A07's Claim 7 asserted

```
UB2_Min_j = min(UB2_j(low), UB2_j(upp))
```

as a computational shortcut for `UB2_Min_j = min_{d in [low,upp]} UB2_j(d)`.
The proof (Remark 5.5.7) has a gap: it only establishes that any critical
point `t* = 1/d*` satisfies `t* < lambda_max(K)` (with
`K = Q^{-1/2} P Q^{-1/2}`, `Q = X'WX`, `P` = coefficient prior precision),
not the claimed `t* < lambda_min(K)`. For anisotropic `K` (coefficient prior
not a Zellner g-prior) this leaves room for a genuine interior minimum below
both endpoints -- confirmed with concrete 2D/3D examples in
`data-raw/ub2_root_finding_prototype.R` (endpoint method overstates the true
minimum by 50-58% there) and, more importantly, with the real compiled
sampler in `data-raw/validate_ub2_rootfinding_fix.R` (see "Validation
results" below). When a proposed/accepted dispersion lands near that
interior point, the code compares against the wrong (too-large) `UB2min_j`
and produces the observed "Sign violation: UB2 < 0" errors.

The only correctness requirement is `UB2_Min_j <= UB2_j(d)` for all `d` in
`[low, upp]` (see Claim 1/3); Claim 7's endpoint formula was never a separate
correctness requirement, just an (invalid, for anisotropic `K`) shortcut for
computing it. Claim 6 (`RSS_Min = min_j RSS_j(low)`), Claims 2, 4, 5, and the
across-face (`PLSD_j`) / across-block (summed `UB2min`) machinery are all
unaffected and require no changes.

## Scope: exactly one file, one function

The fix is fully localized to `bound_ub2_over_dispersion` in
`src/EnvelopeDispersionBuild.cpp` (plus three small, self-contained helper
functions added immediately above it). A byte-level diff confirms
**glmbayes's copy of this file is line-for-line identical** in this function
(glmbayesCore only adds two unrelated diagnostic output fields,
`ub2_at_low`/`ub2_at_upp`, and an `#include <Rmath.h>`) -- so the same patch
is expected to apply to `glmbayes` with no adaptation beyond that.

`minimize_ub2_over_dispersion` (the other function that does endpoint-only
UB2 comparison in the same file) is **dead code in both packages** -- it is
defined but never called anywhere (only `bound_ub2_over_dispersion` is
invoked, from `EnvelopeDispersionBuild_cpp`). It does **not** need to be
touched; leaving it as-is is correct and lower-risk.

## Required changes

### 1. Includes (top of `src/EnvelopeDispersionBuild.cpp`)

Add, next to the existing includes:

```cpp
#include <vector>
#include <algorithm>
#include <cmath>
```

### 2. New helper code -- insert immediately before `bound_ub2_over_dispersion`

Insert the following block verbatim (this is the exact, already-implemented
and validated glmbayesCore code):

```cpp
// ---------------------------------------------------------------------
// Exact (root-finding) minimization of UB2_j(d) over d in [low, upp] for
// anisotropic coefficient priors.
//
// Background (see vignettes/Chapter-A07.Rmd, Remark 5.5.4/5.5.7, and
// data-raw/ub2_root_finding_prototype.R / data-raw/README_ub2_rootfinding_fix.md):
// with t = 1/d, K = Q^{-1/2} P Q^{-1/2}, v_j = Q^{-1/2}*(cbars_j - P*mu - P*beta_hat),
// and w_i = (u_i^T v_j)^2 (u_i = eigenvectors of K, lambda_i its eigenvalues),
//
//   tilde{UB2}_j(t) = (t/2) * (g(t) - Delta),   g(t) = sum_i w_i/(lambda_i+t)^2.
//
// Claim 7 (Chapter A07) assumed the minimum over t always occurs at an
// endpoint, but the underlying proof (Remark 5.5.7) only guarantees any
// critical point t* satisfies t* < lambda_max(K), not t* < lambda_min(K).
// For anisotropic K this allows genuine interior minima, which the
// endpoint-only shortcut misses -- the mechanism behind observed
// "Sign violation: UB2 < 0" errors downstream (the endpoint estimate is too
// large, so an actually-evaluated dispersion near the true interior minimum
// undercuts it). This block finds the true minimum exactly, at negligible
// extra cost: t is always scalar regardless of the coefficient dimension p,
// and K depends only on P/Q (not on the face j), so its eigendecomposition
// is computed once per envelope build and reused across all faces.
// ---------------------------------------------------------------------

namespace ub2_exact_detail {

inline double g_of_t(const arma::vec& lambda, const arma::vec& w, double t) {
  double s = 0.0;
  for (arma::uword i = 0; i < lambda.n_elem; ++i) {
    double m = lambda(i) + t;
    s += w(i) / (m * m);
  }
  return s;
}

inline double hprime_of_t(const arma::vec& lambda, const arma::vec& w, double t) {
  double s = 0.0;
  for (arma::uword i = 0; i < lambda.n_elem; ++i) {
    double m = lambda(i) + t;
    s += w(i) * (lambda(i) - t) / (m * m * m);
  }
  return s;
}

inline double ub2_reduced(const arma::vec& lambda, const arma::vec& w, double Delta, double t) {
  return 0.5 * t * (g_of_t(lambda, w, t) - Delta);
}

// Robust bracketed bisection root finder for f(t) = 0 on [a, b], f(a) and
// f(b) assumed to have opposite signs (or zero). Not the fastest possible
// (Brent's method would converge faster) but simple and reliably correct;
// this runs a handful of times per face, once per envelope build, so
// performance is not a concern.
template <typename F>
inline double bisection_root(F f, double a, double b, double fa, double fb,
                              int max_iter = 100, double tol = 1e-12) {
  for (int it = 0; it < max_iter; ++it) {
    double mid = 0.5 * (a + b);
    double fm = f(mid);
    if (std::abs(fm) < tol || (b - a) < tol * std::max(1.0, std::abs(mid))) return mid;
    if ((fa > 0.0) == (fm > 0.0)) { a = mid; fa = fm; } else { b = mid; fb = fm; }
  }
  return 0.5 * (a + b);
}

struct ExactResult { double ub2_min; double t_star; };

// Finds the exact minimum of tilde{UB2}_j(t) over t in [t_lo, t_hi] by
// bracketing sign changes of h'(t) - Delta on a grid anchored at t_lo, t_hi,
// and min(t_hi, lambda_max(K)) (any critical point must lie strictly below
// lambda_max(K); see Remark 5.5.7), refined near each eigenvalue of K, then
// polishing each bracket via bisection. Always includes t_lo and t_hi among
// the evaluated candidates, so the result is never worse than the old
// endpoint-only estimate.
inline ExactResult ub2_min_exact_1d(const arma::vec& lambda, const arma::vec& w,
                                     double Delta, double t_lo, double t_hi,
                                     int grid_mult = 40) {
  double lam_max = lambda.max();
  double hi_search = std::min(t_hi, lam_max * (1.0 - 1e-9));

  std::vector<double> cands;
  cands.push_back(t_lo);
  cands.push_back(t_hi);

  if (hi_search > t_lo) {
    std::vector<double> anchors;
    anchors.push_back(t_lo);
    anchors.push_back(hi_search);
    for (arma::uword i = 0; i < lambda.n_elem; ++i) {
      if (lambda(i) > t_lo && lambda(i) < hi_search) anchors.push_back(lambda(i));
    }
    std::sort(anchors.begin(), anchors.end());
    anchors.erase(std::unique(anchors.begin(), anchors.end()), anchors.end());

    std::vector<double> grid;
    for (size_t i = 0; i + 1 < anchors.size(); ++i) {
      double lo_i = anchors[i], hi_i = anchors[i + 1];
      if (hi_i <= lo_i) continue;
      for (int k = 0; k < grid_mult; ++k) {
        grid.push_back(lo_i + (hi_i - lo_i) * static_cast<double>(k) / (grid_mult - 1));
      }
    }
    if (grid.size() < 2) {
      for (int k = 0; k < grid_mult; ++k) {
        grid.push_back(t_lo + (hi_search - t_lo) * static_cast<double>(k) / (grid_mult - 1));
      }
    }
    std::sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end()), grid.end());

    std::vector<double> fvals(grid.size());
    for (size_t i = 0; i < grid.size(); ++i) fvals[i] = hprime_of_t(lambda, w, grid[i]) - Delta;

    for (size_t i = 0; i + 1 < grid.size(); ++i) {
      if ((fvals[i] > 0.0) != (fvals[i + 1] > 0.0)) {
        double root = bisection_root(
          [&](double t) { return hprime_of_t(lambda, w, t) - Delta; },
          grid[i], grid[i + 1], fvals[i], fvals[i + 1]
        );
        if (root > t_lo && root < t_hi) cands.push_back(root);
      }
    }
  }

  double best_val = ub2_reduced(lambda, w, Delta, cands[0]);
  double best_t = cands[0];
  for (size_t i = 1; i < cands.size(); ++i) {
    double v = ub2_reduced(lambda, w, Delta, cands[i]);
    if (v < best_val) { best_val = v; best_t = cands[i]; }
  }
  return ExactResult{ best_val, best_t };
}

}  // namespace ub2_exact_detail
```

### 3. Modified function -- `bound_ub2_over_dispersion` body

Find the existing (endpoint-only) body:

```cpp
Rcpp::List bound_ub2_over_dispersion(
    int gs,
    double low,
    double upp,
    const Rcpp::List& cache,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    double rss_min_global
    ) {
  using namespace Rcpp;
  using namespace arma;

  int p = static_cast<int>(cbars.ncol());

  NumericVector disp_min_ub2(gs);
  NumericVector ub2_min(gs);
  NumericVector ub2_at_low(gs);
  NumericVector ub2_at_upp(gs);

  int UB_Min_Method = 2;  // (glmbayes may not have ub2_at_low/ub2_at_upp diagnostics; keep as-is)

  mat base_A  = cache["base_A"];
  vec base_B0 = cache["base_B0"];
  mat Pmat    = cache["Pmat"];
  vec Pmu     = cache["Pmu"];
  Pmat = 0.5 * (Pmat + Pmat.t());

  vec beta_hat;
  mat M_min, M_max;
  double RSS_ML_local = 0.0;

  if (UB_Min_Method == 2) {
    M_min = Rcpp::as<mat>(cache["M_min"]);
    M_max = Rcpp::as<mat>(cache["M_max"]);

    beta_hat      = -solve(base_A, base_B0);
    mat X         = Rcpp::as<mat>(x);
    vec yv        = Rcpp::as<vec>(y);
    vec alphav    = Rcpp::as<vec>(alpha);
    vec wv        = Rcpp::as<vec>(wt);
    vec resid_ml  = yv - X * beta_hat - alphav;
    RSS_ML_local  = as_scalar(resid_ml.t() * (wv % resid_ml));
  }

  for (int j = 0; j < gs; ++j) {
    NumericVector cbars_j(p);
    for (int r = 0; r < p; ++r) cbars_j[r] = cbars(j, r);

    double ub2_low = NA_REAL;
    double ub2_upp = NA_REAL;

    if (UB_Min_Method == 1) {
      ub2_low = UB2(low, cache, cbars_j, y, x, alpha, wt, rss_min_global);
      ub2_upp = UB2(upp, cache, cbars_j, y, x, alpha, wt, rss_min_global);
    } else {
      vec cbar_j_vec(p);
      for (int r = 0; r < p; ++r) cbar_j_vec(r) = cbars(j, r);
      vec b_j = cbar_j_vec - Pmu - Pmat * beta_hat;

      double rss_low_approx = RSS_ML_local + as_scalar(b_j.t() * M_min * b_j);
      double rss_upp_approx = RSS_ML_local + as_scalar(b_j.t() * M_max * b_j);

      ub2_low = (0.5 / low) * (rss_low_approx - rss_min_global);
      ub2_upp = (0.5 / upp) * (rss_upp_approx - rss_min_global);
    }

    ub2_at_low[j] = ub2_low;
    ub2_at_upp[j] = ub2_upp;

    if (ub2_low <= ub2_upp) {
      disp_min_ub2[j] = low;
      ub2_min[j]      = ub2_low;
    } else {
      disp_min_ub2[j] = upp;
      ub2_min[j]      = ub2_upp;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("disp_min_ub2") = disp_min_ub2,
    Rcpp::Named("ub2_min")      = ub2_min,
    Rcpp::Named("ub2_at_low")   = ub2_at_low,
    Rcpp::Named("ub2_at_upp")   = ub2_at_upp
  );
}
```

Replace the local-variable setup and the per-face loop with:

```cpp
Rcpp::List bound_ub2_over_dispersion(
    int gs,
    double low,
    double upp,
    const Rcpp::List& cache,
    const Rcpp::NumericMatrix& cbars,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    double rss_min_global
    ) {
  using namespace Rcpp;
  using namespace arma;

  int p = static_cast<int>(cbars.ncol());

  NumericVector disp_min_ub2(gs);
  NumericVector ub2_min(gs);
  NumericVector ub2_at_low(gs);
  NumericVector ub2_at_upp(gs);

  int UB_Min_Method = 2;

  mat base_A  = cache["base_A"];
  vec base_B0 = cache["base_B0"];
  mat Pmat    = cache["Pmat"];
  vec Pmu     = cache["Pmu"];
  Pmat = 0.5 * (Pmat + Pmat.t());

  vec beta_hat;
  mat M_min, M_max;
  double RSS_ML_local = 0.0;

  // Exact (root-finding) minimization needs K = Q^{-1/2} P Q^{-1/2} and its
  // eigendecomposition. K depends only on P and Q = base_A, not on the face
  // j, so both are computed once here and reused for every face below.
  bool have_K = false;
  vec K_eigval;
  mat K_eigvec;
  mat Qinvhalf;   // Rq^{-1}, where base_A = Rq^T * Rq (arma::chol, upper-tri)
  double Delta = 0.0;
  double t_lo = 1.0 / upp;   // t = 1/d; d=upp -> t_lo (smallest t)
  double t_hi = 1.0 / low;   // d=low -> t_hi (largest t)

  if (UB_Min_Method == 2) {
    M_min = Rcpp::as<mat>(cache["M_min"]);
    M_max = Rcpp::as<mat>(cache["M_max"]);

    beta_hat      = -solve(base_A, base_B0);
    mat X         = Rcpp::as<mat>(x);
    vec yv        = Rcpp::as<vec>(y);
    vec alphav    = Rcpp::as<vec>(alpha);
    vec wv        = Rcpp::as<vec>(wt);
    vec resid_ml  = yv - X * beta_hat - alphav;
    RSS_ML_local  = as_scalar(resid_ml.t() * (wv % resid_ml));

    Delta = rss_min_global - RSS_ML_local;
    if (Delta < 0.0) Delta = 0.0;  // guard against tiny floating-point noise

    mat Rq;
    if (arma::chol(Rq, base_A)) {
      Qinvhalf = arma::inv(arma::trimatu(Rq));
      mat K = Qinvhalf.t() * Pmat * Qinvhalf;
      K = 0.5 * (K + K.t());
      if (arma::eig_sym(K_eigval, K_eigvec, K) && K_eigval.min() > 0.0) {
        have_K = true;
      }
    }
  }

  for (int j = 0; j < gs; ++j) {
    NumericVector cbars_j(p);
    for (int r = 0; r < p; ++r) cbars_j[r] = cbars(j, r);

    double ub2_low = NA_REAL;
    double ub2_upp = NA_REAL;
    double best_ub2, best_disp;
    bool used_exact = false;

    if (UB_Min_Method == 1) {
      // Method 1: original UB2 helper (endpoint-only; retained for
      // compatibility if UB_Min_Method is manually switched back to 1).
      ub2_low = UB2(low, cache, cbars_j, y, x, alpha, wt, rss_min_global);
      ub2_upp = UB2(upp, cache, cbars_j, y, x, alpha, wt, rss_min_global);
    } else {
      // Method 2: RSS-based quadratic form with M_min / M_max
      vec cbar_j_vec(p);
      for (int r = 0; r < p; ++r) cbar_j_vec(r) = cbars(j, r);
      vec b_j = cbar_j_vec - Pmu - Pmat * beta_hat;

      double rss_low_approx = RSS_ML_local + as_scalar(b_j.t() * M_min * b_j);
      double rss_upp_approx = RSS_ML_local + as_scalar(b_j.t() * M_max * b_j);

      ub2_low = (0.5 / low) * (rss_low_approx - rss_min_global);
      ub2_upp = (0.5 / upp) * (rss_upp_approx - rss_min_global);

      if (have_K) {
        // Exact minimum over the whole [low, upp] interval, not just the
        // endpoints -- fixes the Claim 7 gap for anisotropic K. Always
        // evaluates at t_lo/t_hi too, so this is never worse than the
        // endpoint-only estimate above, and reduces to it automatically
        // when K is (numerically) isotropic.
        vec v_j = Qinvhalf.t() * b_j;
        vec w_coords = arma::square(K_eigvec.t() * v_j);
        ub2_exact_detail::ExactResult ex =
          ub2_exact_detail::ub2_min_exact_1d(K_eigval, w_coords, Delta, t_lo, t_hi);
        best_ub2   = ex.ub2_min;
        best_disp  = 1.0 / ex.t_star;
        used_exact = true;
      }
    }

    ub2_at_low[j] = ub2_low;
    ub2_at_upp[j] = ub2_upp;

    if (!used_exact) {
      if (ub2_low <= ub2_upp) {
        best_ub2  = ub2_low;
        best_disp = low;
      } else {
        best_ub2  = ub2_upp;
        best_disp = upp;
      }
    }

    disp_min_ub2[j] = best_disp;
    ub2_min[j]      = best_ub2;
  }

  return Rcpp::List::create(
    Rcpp::Named("disp_min_ub2") = disp_min_ub2,
    Rcpp::Named("ub2_min")      = ub2_min,
    Rcpp::Named("ub2_at_low")   = ub2_at_low,
    Rcpp::Named("ub2_at_upp")   = ub2_at_upp
  );
}
```

If glmbayes's copy does not have the `ub2_at_low`/`ub2_at_upp` diagnostic
fields (glmbayesCore added them separately from this fix), keep whatever
glmbayes already returns for `disp_min_ub2`/`ub2_min` -- only the *values*
assigned to those two need to change, per the diff above; the return-list
shape/field names otherwise stay whatever glmbayes already has.

### 4. Do NOT port the Zellner g-prior guard

An earlier, now-superseded plan (`data-raw/README_glmbayes_g_prior_guard.md`,
kept for historical context but superseded by this document) proposed adding
a `.ing_stop_if_not_g_prior()` R-level guard to reject anisotropic priors as
a stopgap. That guard was implemented, validated, and then **removed** in
glmbayesCore once this root-finding fix was validated directly against it
(see "Validation results" below) -- `R/ing_prior_guard.R` and
`R/simfunction.R`'s call site in glmbayesCore no longer have it. Do not add
it to `glmbayes`; apply the C++ fix above instead, which removes the need
for any such restriction from the start.

## Validation results

Validated in `data-raw/ub2_root_finding_prototype.R` (pure R prototype,
before implementation):
- Reduction formula matches a from-scratch low-level `RSS`/`UB2` replica to
  ~1e-13 (p=2, p=3, real regression data).
- Concrete 2D (`kappa=2000`, 54% gap) and 3D (58% gap) examples with a
  genuine interior minimum, cross-validated against brute-force grid search.
- Condition-number sweep confirms zero gap for `kappa <= 100` and growing
  gap for larger `kappa`.

Validated in `data-raw/validate_ub2_rootfinding_fix.R` (real compiled C++
fix, exercised through the actual `rindepNormalGamma_reg()` sampler with the
temporary g-prior guard disabled): 4 scenarios with deliberately anisotropic,
non-Zellner priors (`K` condition numbers 60,000-800,000, `n` from 6 to 40,
`p` from 2 to 4), 40 repetitions each (160 total), **zero** "Sign violation:
UB2 < 0" errors or diagnostic warnings. The same scenarios reliably produced
these errors before the fix (per the original bug reports that motivated
this work).

## What does NOT change

- `RSS_Min` computation (Claim 6, `bound_rss_over_dispersion`) -- untouched.
- Claims 2, 4, 5 (`test1`, `UB3A`, `UB3B`) -- untouched, don't reference
  `UB2_min_j`.
- The across-face `PLSD_j` mixture-of-truncated-normals proposal -- consumes
  `UB2_min_j` as an opaque per-face constant.
- The across-block summing in `block_rIndepNormalGammaReg.cpp`
  (`UB2min_used += UB2min_j[J_draw[j]]`) -- consumes the (now-corrected)
  per-group `UB2min` values; fixing the one file above fixes both the
  single-group and block/group samplers automatically. No changes needed in
  `rIndepNormalGammaReg.cpp` (confirmed present in both packages, only
  consumes `UB2min` via `UB_list["UB2min"]` in both).
  **Correction:** `block_rIndepNormalGammaReg.cpp` does **not** exist in
  `glmbayes` -- confirmed by directly checking `C:\Rpackages\glmbayes\src`
  (only `rIndepNormalGammaReg.cpp` and `EnvelopeDispersionBuild.cpp` are
  present there). `glmbayes` has no block/group-level sampler at all (no
  `lmbBlock`, no `BlockEnvelopeDispersionBuild`); that feature is exclusive
  to `glmbayesCore` (consumed by `lmebayes`). So there is nothing to port
  for the block sampler in `glmbayes` -- not because it needs no changes,
  but because the file/feature is simply not present there.
- `minimize_ub2_over_dispersion` -- dead code in both packages, not called;
  left as-is.

## Vignette correction needed (glmbayesCore *and* glmbayes, if duplicated)

`vignettes/Chapter-A07.Rmd` has been corrected in glmbayesCore (Remark 5.5.7
and Claim 7); if `glmbayes` ships an analogous copy of this chapter/vignette,
the same textual correction should be applied there. Specifically:

- **Remark 5.5.7**: the proof derives a bound on any critical point `t*` of
  `UB2_j(d)` (in the `t = 1/d` parametrization) and had concluded
  `t* < lambda_min(K)`. The derivation only supports the weaker
  `t* < lambda_max(K)` (the step bounding the critical-point equation used
  `lambda_min(K) <= lambda_i` for all `i`, which gives an upper bound on `t*`
  in terms of `lambda_max`, not `lambda_min` -- `lambda_min(K)` is only a
  valid bound when `K` is isotropic, i.e. `lambda_min = lambda_max`). This
  means the corollary "no critical point can lie in `[lambda_min(K),
  lambda_max(K)]`" is false in general; for anisotropic `K` a critical point
  (and hence an interior local minimum of `UB2_j`) can occur strictly inside
  `(t_lo, t_hi) cap (0, lambda_max(K))`.
- **Claim 7**: currently states `UB2_Min_j = min(UB2_j(low), UB2_j(upp))`
  unconditionally. This should be restated as: this endpoint formula is
  *exact when `K` is isotropic* (Zellner g-prior for the design/weights in
  use); in general, `UB2_Min_j = min_{d in [low,upp]} UB2_j(d)` must be
  computed by minimizing over the interval, e.g. via the root-finding
  recipe in Remark 5.5.4 extended with the corrected bound from 5.5.7
  (scan/bracket/solve `h'(t) = Delta` for `t` in
  `(t_lo, min(t_hi, lambda_max(K)))`, then take the minimum of `tilde{UB2}_j`
  over `{t_lo, t_hi} union {roots}`). The rest of the chapter (Claims 1-6,
  the across-face and across-block machinery) does not depend on Claim 7's
  endpoint shortcut and needs no changes.
