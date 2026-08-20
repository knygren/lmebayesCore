# Numerical checks for inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md.
#
# All checks use a Gaussian-closure setup, where every quantity the note claims
# has an independently known value.
#
# Setup.  Random P11 > 0, random B with 0 < B < P11 (the true S), Pi = P11 - B.
# A deliberately loose floor is built as B_UB = B + t (P11 - B) for t in [0,1),
# giving Pi_LB = P11 - B_UB = (1-t) Pi > 0 -- looseness indexed by one scalar t,
# with t = 0 the tight case.
#
# Checks:
#   1. nu_i = kappa_i^UB / (1 - kappa_i^UB)                        (section 5)
#   2. c' = sqrt(det P11 / det Pi_LB) = prod (1-kappa^UB)^(-1/2)   (section 4)
#   3. c' = M_Q(1) = E_Q[exp(Psi)] under closure                   (section 4)
#   4. sandwich 1 <= c <= c'                                       (section 4)
#   5. M* = c exp(d - r_d/2) with r_d >= 2d; M* = c = 1 at t = 0   (section 2)
#   6. M* is the sup of pi/q2 over Ctilde^c (boundary attainment)  (section 2)
#   7. the section 5 bound dominates the exact escape probability
#
# Weighted-chisq tails are evaluated by Monte Carlo with COMMON RANDOM NUMBERS
# across the bound and the exact probability, so that the tight case (t = 0)
# agrees identically rather than up to simulation error.
#
# Scratch check, not part of the package.

set.seed(20260816)

sqrtm <- function(M) {
  e <- eigen(M, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(M)) %*% t(e$vectors)
}
isqrtm <- function(M) solve(sqrtm(M))
gen_eig <- function(X, Y) sort(Re(eigen(solve(Y, X), only.values = TRUE)$values))

qmax <- 5L
nmc  <- 1e6L
Z2   <- matrix(rnorm(nmc * qmax), ncol = qmax)^2   # shared across all cases

wtail <- function(x, lambda) mean(Z2[, seq_along(lambda), drop = FALSE] %*%
                                    lambda > x)

make_case <- function(q, t) {
  U   <- qr.Q(qr(matrix(rnorm(q * q), q)))
  P11 <- U %*% diag(sort(runif(q, 1, 5)), q) %*% t(U)
  V   <- qr.Q(qr(matrix(rnorm(q * q), q)))
  D   <- V %*% diag(runif(q, 0.1, 0.9), q) %*% t(V)   # 0 < D < I
  R   <- sqrtm(P11)
  B   <- (function(M) (M + t(M)) / 2)(R %*% D %*% R)
  Pi  <- P11 - B
  list(q = q, t = t, P11 = P11, B = B, Pi = Pi,
       B_UB = B + t * Pi, Pi_LB = P11 - (B + t * Pi))
}

cases <- c(
  lapply(1:60, function(i) make_case(sample(2:qmax, 1), runif(1, 0.05, 0.9))),
  lapply(1:10, function(i) make_case(sample(2:qmax, 1), 0))
)

d <- 3

res <- do.call(rbind, lapply(cases, function(cs) {
  kUB <- gen_eig(cs$B_UB, cs$P11)
  nu  <- gen_eig(cs$B_UB, cs$Pi_LB)
  mu  <- gen_eig(cs$B,    cs$Pi)

  cp_det  <- sqrt(det(cs$P11) / det(cs$Pi_LB))
  cp_prod <- prod((1 - kUB)^(-1 / 2))
  cc      <- sqrt(det(cs$Pi) / det(cs$Pi_LB))     # closure value of c

  Bi  <- isqrtm(cs$B)
  lam <- min(eigen(Bi %*% cs$B_UB %*% Bi, symmetric = TRUE,
                   only.values = TRUE)$values)
  r_d <- 2 * d * lam
  Mst <- cc * exp(d - r_d / 2)

  bound <- Mst * wtail(2 * d, nu)
  exact <- wtail(2 * d, mu)

  data.frame(q = cs$q, t = cs$t,
             err_nu = max(abs(nu - kUB / (1 - kUB))),
             err_cp = abs(cp_det - cp_prod) / cp_prod,
             sandwich_ok = (cc >= 1 - 1e-10) && (cc <= cp_det + 1e-10),
             c = cc, cprime = cp_det, Mstar = Mst,
             r_side_ok = r_d >= 2 * d - 1e-9,
             tight_ok = if (cs$t == 0) abs(Mst - 1) < 1e-8 else NA,
             bound = bound, exact = exact,
             dominates = bound >= exact - 1e-12)
}))

cat("=== Checks 1, 2, 4, 5, 7 over", nrow(res), "closure cases ===\n")
cat("max |nu - kUB/(1-kUB)|:                ", format(max(res$err_nu), digits = 3), "\n")
cat("max rel err in c' (det vs product):    ", format(max(res$err_cp), digits = 3), "\n")
cat("sandwich 1 <= c <= c' holds:           ", all(res$sandwich_ok), "\n")
cat("r_d >= 2d holds:                       ", all(res$r_side_ok), "\n")
cat("M* = c = 1 in all tight (t = 0) cases: ",
    all(res$tight_ok[!is.na(res$tight_ok)]), "\n")
## At t = 0 the two weight vectors are the same pencil computed two ways, so
## they differ only by eigensolver round-off; the MC indicator can then flip for
## a handful of the nmc draws.  Report the magnitude rather than exact equality.
cat("tight cases: max |bound - exact|:      ",
    format(max(abs(res$bound[res$t == 0] - res$exact[res$t == 0])),
           digits = 3), "\n")
cat("bound dominates exact:                 ", all(res$dominates),
    " (violations:", sum(!res$dominates), ")\n")
cat("bound/exact ratio, loose cases:        ",
    format(range(res$bound[res$t > 0] / res$exact[res$t > 0]), digits = 4), "\n")
cat("M*/c ratio, loose cases:               ",
    format(range(res$Mstar[res$t > 0] / res$c[res$t > 0]), digits = 3), "\n")
cat("c'/c ratio, loose cases:               ",
    format(range(res$cprime[res$t > 0] / res$c[res$t > 0]), digits = 3), "\n")

cat("\n=== Check 3: c' equals M_Q(1) = E_Q[exp(Psi)] under closure ===\n")
cs  <- make_case(3, 0)
n3  <- 2e6
X   <- matrix(rnorm(n3 * 3), ncol = 3) %*% chol(solve(cs$P11))
Psi <- 0.5 * rowSums((X %*% cs$B) * X)
mq1 <- mean(exp(Psi))
cpr <- sqrt(det(cs$P11) / det(cs$Pi_LB))
cat(sprintf("MC E_Q[exp(Psi)] = %.6f   c' = %.6f   rel err = %.2e\n",
            mq1, cpr, abs(mq1 - cpr) / cpr))

## Check 9 (Theorem 2 trade-off for logit) lives in its own file, since it needs
## a simulated data set and full quadrature rather than a closed-form identity:
## see data-raw/_ex_logit_majorization_floor.R, written up in sections 3A.6-3A.7.

cat("\n=== Check 8: S is the Omega-weighted P11_RE (section 6) ===\n")
## S = sum_j H_j' Pb^{1/2} Omega_j Pb^{1/2} H_j with
## Omega_j = (I + Pb^{-1/2} G_j Pb^{-1/2})^{-1}, and
## kappa_0 = lmax(P11^{-1} S) <= omega_max * lmax(P11^{-1} P11_RE).
err_S <- err_k <- numeric(0)
for (rep in 1:200) {
  q  <- sample(2:4, 1); pre <- sample(1:3, 1); J <- sample(2:6, 1)
  Pb <- (function(M) crossprod(M) + diag(pre)) (matrix(rnorm(pre * pre), pre))
  Rb <- sqrtm(Pb); Rbi <- solve(Rb)
  H  <- lapply(1:J, function(j) matrix(rnorm(pre * q), pre, q))
  G  <- lapply(1:J, function(j) crossprod(matrix(rnorm(4 * pre), 4, pre)))
  Lg <- (function(M) crossprod(M) + diag(q) * 0.3)(matrix(rnorm(q * q), q))

  P11_RE <- Reduce(`+`, lapply(1:J, function(j) t(H[[j]]) %*% Pb %*% H[[j]]))
  P11    <- Lg + P11_RE
  S_dir  <- Reduce(`+`, lapply(1:J, function(j)
    t(H[[j]]) %*% Pb %*% solve(Pb + G[[j]]) %*% Pb %*% H[[j]]))
  Om     <- lapply(1:J, function(j) solve(diag(pre) + Rbi %*% G[[j]] %*% Rbi))
  S_om   <- Reduce(`+`, lapply(1:J, function(j)
    t(H[[j]]) %*% Rb %*% Om[[j]] %*% Rb %*% H[[j]]))

  err_S <- c(err_S, max(abs(S_dir - S_om)) / max(abs(S_dir)))

  om_max <- max(sapply(Om, function(O)
    max(eigen(O, symmetric = TRUE, only.values = TRUE)$values)))
  k0     <- max(gen_eig(S_dir, P11))
  k_bnd  <- om_max * max(gen_eig(P11_RE, P11))
  err_k  <- c(err_k, k0 - k_bnd)
}
cat("max rel err, S direct vs Omega form:   ", format(max(err_S), digits = 3), "\n")
cat("kappa_0 <= omega_max * lmax(P11^-1 P11_RE) holds:",
    all(err_k <= 1e-10), " (worst slack:",
    format(max(err_k), digits = 3), ")\n")

cat("\n=== Check 6: boundary attainment of M* (q = 2 grid scan) ===\n")
for (tt in c(0.1, 0.4, 0.8)) {
  cs <- make_case(2, tt)
  cc <- sqrt(det(cs$Pi) / det(cs$Pi_LB))
  g  <- seq(-40, 40, length.out = 2001)
  G  <- as.matrix(expand.grid(g, g))
  ext <- rowSums((G %*% cs$B) * G) > 2 * d
  Rv  <- -0.5 * rowSums((G %*% (cs$B_UB - cs$B)) * G)
  obs <- cc * max(exp(Rv[ext]))
  Bi  <- isqrtm(cs$B)
  lam <- min(eigen(Bi %*% cs$B_UB %*% Bi, symmetric = TRUE,
                   only.values = TRUE)$values)
  cat(sprintf("t=%.1f  grid sup=%.8f  closed form=%.8f  rel err=%.2e\n",
              tt, obs, cc * exp(d - d * lam),
              abs(obs - cc * exp(d - d * lam)) / (cc * exp(d - d * lam))))
}
