# SUPERSEDED -- do not port this guard to `glmbayes`

**Status: superseded.** This guard was implemented and tested in
**glmbayesCore**, but has since been **removed** there: the root cause it
worked around (Chapter A07 Claim 7's endpoint-only `UB2_Min_j` being invalid
for anisotropic coefficient priors) is now fixed directly via exact
root-finding in `src/EnvelopeDispersionBuild.cpp::bound_ub2_over_dispersion`.
See `data-raw/README_ub2_rootfinding_fix.md` for the fix and the porting
instructions for `glmbayes` -- **apply that fix instead of this guard**. This
file is kept only for historical context (the guard was a real, working
stopgap for one release of glmbayesCore) and should not be used as a porting
guide anymore.

---

Original content below (historical; do not act on it).

Implemented and tested in **glmbayesCore** (`R/ing_prior_guard.R`,
`R/simfunction.R`). Two changes to port into **glmbayes**'s own
(independent) copy of `rindepNormalGamma_reg()` in `R/simfunction.R`.

## 1. New function — add to `glmbayes/R/simfunction.R` (or a new file)

```r
#' Stop when an ING coefficient prior is not (numerically) a Zellner g-prior
#'
#' TEMPORARY safeguard. See glmbayesCore R/ing_prior_guard.R for full
#' rationale and data-raw/README_g_prior_safeguard.md for background.
#' @noRd
.ing_stop_if_not_g_prior <- function(P, x, wt, rel_tol = 1e-6, prefix = NULL) {
  p <- ncol(x)
  if (p < 2L) {
    return(invisible(1))
  }

  Q <- crossprod(x * sqrt(wt))
  Q <- 0.5 * (Q + t(Q))

  Rq <- tryCatch(
    chol(Q),
    error = function(e) {
      stop(
        prefix,
        "cannot verify the Zellner g-prior restriction: t(x) %*% diag(wt) %*% x ",
        "is not positive definite (", conditionMessage(e), ").",
        call. = FALSE
      )
    }
  )
  A <- backsolve(Rq, diag(p))
  K <- crossprod(A, P %*% A)
  K <- 0.5 * (K + t(K))

  ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  lambda_min <- min(ev)
  lambda_max <- max(ev)

  if (!is.finite(lambda_min) || !is.finite(lambda_max) || lambda_min <= 0) {
    stop(
      prefix,
      "cannot verify the Zellner g-prior restriction: K = Q^{-1/2} P Q^{-1/2} ",
      "is not positive definite.",
      call. = FALSE
    )
  }

  ratio <- lambda_max / lambda_min

  if (ratio - 1 > rel_tol) {
    stop(
      prefix,
      "dIndependent_Normal_Gamma currently requires the coefficient prior ",
      "'Sigma' to be (numerically) a Zellner g-prior for THIS call's design ",
      "and weights, i.e. Sigma proportional to ",
      "solve(t(x) %*% diag(wt) %*% x). The supplied Sigma is not: ",
      "K = Q^{-1/2} P Q^{-1/2} has eigenvalue ratio lambda_max/lambda_min = ",
      signif(ratio, 4), " (expected ~1 within ", signif(rel_tol, 2), "). ",
      "This restriction is temporary; the accept-reject envelope's ",
      "dispersion-bound derivation (Chapter A07 Claim 7) has a gap for ",
      "anisotropic priors and can otherwise trigger spurious 'UB2 < 0' ",
      "sign-violation errors. Use Prior_Setup()'s default (uncustomized) ",
      "Sigma.",
      call. = FALSE
    )
  }

  invisible(ratio)
}
```

## 2. Modified function — `rindepNormalGamma_reg()` in `glmbayes/R/simfunction.R`

Insert one guard call right after the existing `P` SPD check (line ~1242
in the current file):

```1237:1243:c:\Rpackages\glmbayes\R\simfunction.R
  stopifnot(isSymmetric(P))
  
  tol <- 1e-6
  ev  <- eigen(P, symmetric = TRUE)$values
  stopifnot(all(ev >= -tol * abs(ev[1L])))
  
  # dispersion must be numeric scalar or NULL
```

becomes:

```r
  stopifnot(isSymmetric(P))
  
  tol <- 1e-6
  ev  <- eigen(P, symmetric = TRUE)$values
  stopifnot(all(ev >= -tol * abs(ev[1L])))
  
  ## TEMPORARY: require a (numerically) Zellner g-prior for the coefficient
  ## covariance (see .ing_stop_if_not_g_prior() above).
  .ing_stop_if_not_g_prior(P = P, x = x, wt = wt)
  
  # dispersion must be numeric scalar or NULL
```

## Known test impact

`tests/testthat/test-lmb-non-zellner.R` deliberately uses
`Sigma_non_zellner <- 0.001 * diag(diag(ps$Sigma))` and will start failing
(erroring) once this guard is in place, since dropping the off-diagonal
terms of `ps$Sigma` breaks proportionality to `Q^{-1}` (confirmed in
glmbayesCore with an analogous non-Zellner case, rejected for eigenvalue
ratio > 1). That test needs to be updated to `expect_error(...)` (or
removed/retired) as part of this port.
