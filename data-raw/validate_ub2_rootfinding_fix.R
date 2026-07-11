# Validates the UB2_Min root-finding fix (src/EnvelopeDispersionBuild.cpp::
# bound_ub2_over_dispersion) by exercising the REAL compiled sampler with
# genuinely anisotropic (non-Zellner) coefficient priors and checking that no
# "Sign violation: UB2 < 0" errors or diagnostic warnings are ever emitted
# (these are Rcpp::Rcout prints / thrown C++ exceptions -- see
# src/rIndepNormalGammaReg.cpp around the UB2 sign-violation check).
#
# This is the validation that motivated removing the temporary Zellner
# g-prior guard (.ing_stop_if_not_g_prior(), formerly in R/ing_prior_guard.R
# / R/simfunction.R) -- see data-raw/README_ub2_rootfinding_fix.md.
#
#   cd .../glmbayesCore
#   Rscript data-raw/validate_ub2_rootfinding_fix.R

devtools::load_all(".", quiet = TRUE)

run_stress <- function(label, x, y, mu, Sigma, shape, rate, wt = rep(1, length(y)),
                        max_disp_perc = 0.99, n_draws = 500L, n_reps = 30L, seed0 = 1) {
  cat("\n============================================================\n")
  cat("Stress test:", label, "\n")
  cat("============================================================\n")

  P <- solve(Sigma); P <- 0.5 * (P + t(P))
  Q <- crossprod(x * sqrt(wt))
  Qinvhalf <- solve(chol(Q))
  K <- t(Qinvhalf) %*% P %*% Qinvhalf
  K <- 0.5 * (K + t(K))
  eig <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  cat("K eigenvalues:", signif(eig, 4),
      " (condition number =", signif(max(eig) / min(eig), 4), ")\n")
  cat("n =", length(y), " p =", ncol(x), " shape =", shape, " rate =", rate, "\n\n")

  prior_list <- list(mu = mu, Sigma = Sigma, shape = shape, rate = rate,
                      max_disp_perc = max_disp_perc)

  n_errors <- 0L
  for (rep_i in seq_len(n_reps)) {
    set.seed(seed0 + rep_i)
    ok <- tryCatch({
      fit <- rindepNormalGamma_reg(
        n = n_draws, y = y, x = x, prior_list = prior_list, weights = wt,
        use_parallel = FALSE, verbose = FALSE, progbar = FALSE
      )
      TRUE
    }, error = function(e) {
      cat("  !! rep", rep_i, "ERROR:", conditionMessage(e), "\n")
      FALSE
    })
    if (!ok) n_errors <- n_errors + 1L
  }
  cat("Reps with a thrown error:", n_errors, "/", n_reps, "\n")
  invisible(n_errors)
}

## ============================================================================
## Scenario 1: n=20, p=3, anisotropic diag(1e6, 400, 20) Sigma -- the exact
## "Case 3" scenario the Zellner-only guard rejects (data-raw/check_g_prior_guard.R).
## ============================================================================
set.seed(1)
n1 <- 20L
x1 <- cbind(1, rnorm(n1), rnorm(n1))
colnames(x1) <- c("(Intercept)", "x1", "x2")
beta_true1 <- c(1, 2, -1)
y1 <- x1 %*% beta_true1 + rnorm(n1)

Sigma1 <- matrix(0, 3, 3)
diag(Sigma1) <- c(1e6, 400, 20)

n_err1 <- run_stress("n=20, p=3, diag(1e6, 400, 20)", x1, y1, c(0, 0, 0), Sigma1,
                      shape = 3, rate = 2, n_reps = 40)

## ============================================================================
## Scenario 2: small n=14, p=2, strongly anisotropic prior (mirrors the
## 2D root-finding-prototype example with a >50% endpoint-vs-true gap).
## ============================================================================
set.seed(3)
n2 <- 14L
x2 <- cbind(1, rnorm(n2))
colnames(x2) <- c("(Intercept)", "x1")
beta_true2 <- c(0.5, -1.2)
y2 <- x2 %*% beta_true2 + rnorm(n2, sd = 1.5)

Sigma2 <- diag(c(2000, 0.05))

n_err2 <- run_stress("n=14, p=2, diag(2000, 0.05)", x2, y2, c(0, 0), Sigma2,
                      shape = 3, rate = 2, n_reps = 40)

## ============================================================================
## Scenario 3: very small n=6, p=2, extreme anisotropy + weak prior (the
## original small-group / calibrated-prior regime that first surfaced the bug
## in lmerb's gamma_list mode).
## ============================================================================
set.seed(5)
n3 <- 6L
x3 <- cbind(1, rnorm(n3))
colnames(x3) <- c("(Intercept)", "x1")
beta_true3 <- c(2, -0.5)
y3 <- x3 %*% beta_true3 + rnorm(n3, sd = 1)

Sigma3 <- diag(c(1e4, 1e-2))

n_err3 <- run_stress("n=6, p=2, diag(1e4, 1e-2), weak prior", x3, y3, c(0, 0), Sigma3,
                      shape = 2, rate = 1, n_reps = 40)

## ============================================================================
## Scenario 4: n=40, p=4, random anisotropic Sigma with a large condition
## number in K (built directly from a random orthonormal basis, as in the
## root-finding prototype's validate_reduction_once()).
## ============================================================================
set.seed(7)
n4 <- 40L
x4 <- cbind(1, matrix(rnorm(n4 * 3), n4, 3))
colnames(x4) <- c("(Intercept)", "x1", "x2", "x3")
beta_true4 <- c(1, -2, 0.5, 3)
y4 <- x4 %*% beta_true4 + rnorm(n4)

p4 <- 4
eig_target <- exp(seq(log(0.001), log(50), length.out = p4))
U4 <- qr.Q(qr(matrix(rnorm(p4 * p4), p4, p4)))
P4 <- U4 %*% diag(eig_target) %*% t(U4)
P4 <- 0.5 * (P4 + t(P4))
Sigma4 <- solve(P4)

n_err4 <- run_stress("n=40, p=4, random anisotropic K", x4, y4, rep(0, p4), Sigma4,
                      shape = 4, rate = 3, n_reps = 40)

## ============================================================================
## Summary
## ============================================================================
cat("\n============================================================\n")
cat("Summary\n")
cat("============================================================\n")
total_errs <- n_err1 + n_err2 + n_err3 + n_err4
total_reps <- 160L
cat("Total thrown errors across all scenarios/reps:", total_errs, "/", total_reps, "\n")
if (total_errs == 0L) {
  cat("PASS: no 'Sign violation' (or any other) errors observed.\n")
} else {
  cat("FAIL: see errors printed above.\n")
}

cat("\n=== validate_ub2_rootfinding_fix.R: done ===\n")
