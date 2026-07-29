## TEMPORARY / SCRATCH -- investigation only, not package code.
##
## Reproducible numeric check backing inst/omega-ing-marginal-multivariate-t.md:
## verifies that Section 16.2's closed-form H_j(beta_j) (i) matches a
## finite-difference Hessian of -log p(y_j|beta_j), (ii) matches the general
## multivariate-t Hessian formula from inst/multivariate-t-log-concavity.md
## once nu_j = 2*a0+n-p and Sigma_j = (2*r0+RSS_ols)/nu_j*(D'D)^{-1} are
## substituted, and (iii) that varying the Gamma prior's SHAPE (a0) alone
## rescales H_j's eigenvalues but never flips their sign (consistent with a0
## dropping out of the q_j <= 2*r0+RSS_ols log-concavity boundary exactly).
##
## Also verifies (Sections 7-8 of the same README): (iv) the exact F-to-Beta
## pivot alpha_j = P(F_{p,nu} > nu/p) = P(Beta(p/2, nu/2) > 1/2), and its
## monotonic decrease in nu; (v) that sweeping pwt_measurement moves
## shape_ING,j/rate,j/nu_j/K_j together per dGamma_list()'s own calibration
## formulas (inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md Part I), with alpha_j
## still decreasing monotonically throughout; and (vi) that the generalized-t
## KERNEL converges pointwise to the Gaussian kernel exp(-d/2) as nu -> Inf
## (the "sharper prior -> more Gaussian-like, globally log-concave" limit).

set.seed(1)
n <- 7; p <- 2
D <- matrix(rnorm(n*p), n, p)
beta_true <- c(1, -0.5)
y <- D %*% beta_true + rnorm(n, sd = 0.7)
a0 <- 2.3; r0 <- 1.7

DtD <- crossprod(D)
beta_ols <- as.vector(solve(DtD, crossprod(D, y)))
RSS_ols <- sum((y - D %*% beta_ols)^2)

neglogp <- function(beta) {
  e <- y - D %*% beta
  ete <- sum(e^2)
  (a0 + n / 2) * log(r0 + 0.5 * ete)
}

H_closed <- function(beta) {
  e <- as.vector(y - D %*% beta)
  ete <- sum(e^2)
  Dte <- as.vector(crossprod(D, e))
  Omega_eff <- (a0 + n / 2) / (r0 + 0.5 * ete)
  Omega_eff * DtD - (Omega_eff^2 / (a0 + n / 2)) * outer(Dte, Dte)
}

numHess <- function(f, x, h = 1e-5) {
  p <- length(x)
  H <- matrix(0, p, p)
  for (i in 1:p) for (j in 1:p) {
    e_i <- rep(0, p); e_i[i] <- h
    e_j <- rep(0, p); e_j[j] <- h
    H[i, j] <- (f(x + e_i + e_j) - f(x + e_i - e_j) - f(x - e_i + e_j) + f(x - e_i - e_j)) / (4 * h^2)
  }
  H
}

beta_test <- beta_ols + c(0.3, -0.2)
Hn <- numHess(neglogp, beta_test)
Hc <- H_closed(beta_test)
cat("max abs diff closed vs numeric Hessian:", max(abs(Hn - Hc)), "\n")

nu <- 2 * a0 + n - p
K  <- 2 * r0 + RSS_ols
A  <- DtD * nu / K
z <- beta_test - beta_ols
d <- as.numeric(t(z) %*% A %*% z)
s <- 1 + d / nu
Aq <- A %*% z
Ht <- -((nu + p) / (nu * s)) * A + (2 * (nu + p) / (nu^2 * s^2)) * (Aq %*% t(Aq))
cat("max abs diff generalized-t (-Ht) vs closed H_j:", max(abs((-Ht) - Hc)), "\n")

q_j <- as.numeric(t(z) %*% DtD %*% z)
threshold <- 2 * r0 + RSS_ols
cat("q_j:", q_j, " threshold:", threshold, " (should be independent of a0)\n")
for (a0_try in c(0.01, 1, 5, 50)) {
  ete_test <- sum((y - D %*% beta_test)^2)
  Omega_eff_try <- (a0_try + n / 2) / (r0 + 0.5 * ete_test)
  e <- as.vector(y - D %*% beta_test)
  Dte <- as.vector(crossprod(D, e))
  H_try <- Omega_eff_try * DtD - (Omega_eff_try^2 / (a0_try + n / 2)) * outer(Dte, Dte)
  cat("a0=", a0_try, " min eigenvalue of H_j:",
      min(eigen(H_try, symmetric = TRUE, only.values = TRUE)$values), "\n")
}

## ---------------------------------------------------------------------------
## Section 7 check: alpha_j = P(F_{p,nu} > nu/p) == P(Beta(p/2,nu/2) > 1/2)
## exactly, and is monotonically decreasing in nu (fixed p).
## ---------------------------------------------------------------------------
cat("\n--- Section 7: F-to-Beta pivot + monotonicity in nu ---\n")
p_re <- 3
nu_grid <- c(p_re, 2 * p_re, 5 * p_re, 20 * p_re, 200 * p_re)
alpha_F    <- sapply(nu_grid, function(nu) stats::pf(nu / p_re, p_re, nu, lower.tail = FALSE))
alpha_Beta <- sapply(nu_grid, function(nu) stats::pbeta(0.5, p_re / 2, nu / 2, lower.tail = FALSE))
cat("max abs diff alpha_j via pf() vs pbeta() pivot:", max(abs(alpha_F - alpha_Beta)), "\n")
cat("nu:", nu_grid, "\n")
cat("alpha_j (should be strictly decreasing, = 0.5 at nu = p_re):", alpha_F, "\n")
cat("strictly decreasing:", all(diff(alpha_F) < 0), "\n")
cat("alpha_j at nu = p_re equals 0.5:", isTRUE(all.equal(alpha_F[1], 0.5)), "\n")

## ---------------------------------------------------------------------------
## Section 8 check: sweep pwt_measurement (w_j), rebuild shape_ING,j/rate,j
## from dGamma_list()'s own calibration formulas (DGAMMA_LIST_MARGINAL_AND_
## BOUNDS.md Part I), and confirm a0,r0,nu,K,K/nu increase while alpha_j
## decreases monotonically, and that r0 = sigma2_hat*(a0-1) exactly throughout.
## ---------------------------------------------------------------------------
cat("\n--- Section 8: sweeping pwt_measurement ---\n")
n_j <- 25; p_re2 <- 2
RSS_ols_j   <- 18.4                 # arbitrary fixed group-level quantities
penalty     <- 6.1                  # (bhat_j - mu_j)' Mi_j (bhat_j - mu_j), >= 0
S_marg_j    <- RSS_ols_j + penalty
sigma2_hat_j <- S_marg_j / (n_j - p_re2)

w_grid <- c(0, 0.2, 0.5, 0.8, 0.95, 0.99)
res <- do.call(rbind, lapply(w_grid, function(w) {
  n_prior_j <- if (w == 0) 0 else w / (1 - w) * n_j
  a0j <- (n_prior_j + 1) / 2 + p_re2 / 2
  r0j <- 0.5 * S_marg_j * (n_prior_j + p_re2 - 1) / (n_j - p_re2)
  nuj <- 2 * a0j + n_j - p_re2
  Kj  <- 2 * r0j + RSS_ols_j
  alphaj <- stats::pf(nuj / p_re2, p_re2, nuj, lower.tail = FALSE)
  data.frame(w = w, n_prior = n_prior_j, a0 = a0j, r0 = r0j, nu = nuj,
             K = Kj, K_over_nu = Kj / nuj, alpha = alphaj,
             r0_identity_check = r0j - sigma2_hat_j * (a0j - 1))
}))
print(res, row.names = FALSE, digits = 5)
cat("a0 strictly increasing:", all(diff(res$a0) > 0), "\n")
cat("r0 strictly increasing:", all(diff(res$r0) > 0), "\n")
cat("nu strictly increasing:", all(diff(res$nu) > 0), "\n")
cat("K strictly increasing:", all(diff(res$K) > 0), "\n")
cat("K/nu strictly increasing:", all(diff(res$K_over_nu) > 0), "\n")
cat("alpha strictly decreasing:", all(diff(res$alpha) < 0), "\n")
cat("max abs r0 identity residual (should be ~0):", max(abs(res$r0_identity_check)), "\n")
cat("K/nu approaching sigma2_hat_j =", sigma2_hat_j, "from below as w -> 1\n")

## ---------------------------------------------------------------------------
## Section 8 (cont'd): as nu -> Inf (pwt_measurement -> 1), the generalized-t
## KERNEL [1+d/nu]^{-(nu+p)/2} converges pointwise to the Gaussian kernel
## exp(-d/2) -- the "deeper mechanism" behind alpha_j -> 0 (a Gaussian is
## globally log-concave; inst/multivariate-t-log-concavity.md's own
## "Gaussian limit" remark).
## ---------------------------------------------------------------------------
cat("\n--- Section 8 (cont'd): t-kernel -> Gaussian-kernel as nu -> Inf ---\n")
d_test <- 5
nu_seq <- c(5, 20, 100, 1000, 1e5)
kernel_t     <- sapply(nu_seq, function(nu) (1 + d_test / nu)^(-(nu + p_re2) / 2))
kernel_gauss <- exp(-d_test / 2)
cat("d =", d_test, " p_re =", p_re2, "\n")
cat("nu:            ", nu_seq, "\n")
cat("t-kernel/gauss-kernel:", kernel_t / kernel_gauss, " (-> 1 as nu grows)\n")
