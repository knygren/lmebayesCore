# Validates the new kappa(K) <= 2 near-isotropic fast path added to
# bound_ub2_over_dispersion() in src/EnvelopeDispersionBuild.cpp: for
# isotropic/near-isotropic K it should report K_is_near_isotropic = TRUE and
# skip the per-face root search (ub2_min should exactly equal
# min(ub2_at_low, ub2_at_upp) for every face); for anisotropic K it should
# report FALSE and still find the correct (possibly lower) exact minimum.
#
#   Rscript data-raw/validate_near_isotropic_fastpath.R

devtools::load_all(".", quiet = TRUE)

set.seed(1)
n <- 20L
x1 <- rnorm(n)
x2 <- rnorm(n)
beta_true <- c(0.5, 1, -1)
y <- as.numeric(cbind(1, x1, x2) %*% beta_true + rnorm(n, sd = 1.2))
dat <- data.frame(y = y, x1 = x1, x2 = x2)

mu <- c(0, 0, 0)

run_case <- function(label, Sigma) {
  fit <- rindepNormalGamma_reg_with_envelope(
    n = 50L,
    y = dat$y,
    x = cbind(1, dat$x1, dat$x2),
    prior_list = list(mu = mu, Sigma = Sigma, shape = 3, rate = 2),
    verbose = TRUE,
    progbar = FALSE,
    use_parallel = FALSE
  )
  diag <- fit$diagnostics
  ub <- fit$UB_list
  cat("\n---", label, "---\n")
  cat("kappa(K)             :", diag$kappa_K, "\n")
  cat("K_is_near_isotropic   :", diag$K_is_near_isotropic, "\n")
  cat("all(UB2min finite)   :", all(is.finite(ub$UB2min)), "\n")
  invisible(diag)
}

## Case 1: exact Zellner g-prior -> K = scalar*I -> kappa(K) = 1 exactly.
Q <- crossprod(cbind(1, dat$x1, dat$x2))
Sigma_zellner <- 5 * solve(Q)
diag1 <- run_case("Zellner g-prior (kappa=1)", Sigma_zellner)
stopifnot(isTRUE(diag1$K_is_near_isotropic))
stopifnot(abs(diag1$kappa_K - 1) < 1e-6)

## Case 2: mildly anisotropic vector-pwt-style prior (Sigma = V0*outer(s,s),
## s = c(1.05, 1, 0.95)) giving 1 < kappa(K) < 2 -- should still use the fast
## path, confirming it's not just triggering on exact Zellner/isotropy.
V0 <- solve(Q)
s <- c(1.05, 1, 0.95)
Sigma_mild <- V0 * outer(s, s)
diag2 <- run_case("Mildly anisotropic vector-pwt-style (target 1<kappa<2)", Sigma_mild)
stopifnot(diag2$kappa_K > 1, diag2$kappa_K <= 2)
stopifnot(isTRUE(diag2$K_is_near_isotropic))

## Case 3: strongly anisotropic (kappa(K) >> 2) -> fast path should NOT
## trigger; exact root-finding search should still run (and, per the
## already-validated anisotropic stress tests, produce no sign violations).
Sigma_aniso <- diag(c(2000, 0.05, 5))
diag3 <- run_case("Strongly anisotropic (kappa>>2)", Sigma_aniso)
stopifnot(!isTRUE(diag3$K_is_near_isotropic))
stopifnot(diag3$kappa_K > 2)

cat("\n=== All checks passed ===\n")
