## Scratch check for inst/RESTRICTED_GIBBS_MINORIZATION_TV.md, Lemma 2A.3-F.
## Verifies, on a small non-Gaussian (logit) GLMM by brute-force quadrature:
##   (1) gamma* = argmin Phi = fixed point of M
##   (2) D(gamma) - D(gamma*) = Phibar(gamma) - 0.5 * P11 * (gamma - gamma*)^2
##   (3) Psi := D(gamma*) - D(gamma) satisfies 0 <= Psi <= 0.5 * P11 * (gamma-gamma*)^2
##   (4) Psi is convex and epsilon peaks at gamma*
## Not part of the package; delete freely.

set.seed(1)

## ---- model: q = 1, J = 2 groups, W_j = 1, beta_j | gamma ~ N(gamma, psi) ----
J     <- 2
psi   <- 0.7
Pb    <- 1 / psi
lamg  <- 0.35          # population prior precision Lambda_gamma
mu0   <- 0.2
n_tr  <- c(6, 4)       # binomial trials per group
y_su  <- c(5, 1)       # successes per group

P11re <- J * Pb
P11   <- lamg + P11re

loglik <- function(b, j) y_su[j] * b - n_tr[j] * log1p(exp(b))

## quadrature grid in beta-space
bg <- seq(-14, 14, length.out = 1201)
db <- bg[2] - bg[1]

## unnormalised pi(beta_j | gamma, y) for one group
lw_group <- function(gam, j) loglik(bg, j) - 0.5 * Pb * (bg - gam)^2

## m(beta) = P11^{-1} (lamg*mu0 + sum_j Pb * beta_j)
m_of <- function(b1, b2) (lamg * mu0 + Pb * (b1 + b2)) / P11

## ---- Phi(gamma) via the log-partition A_j ----
A_j <- function(theta, j) {
  lw <- loglik(bg, j) - 0.5 * Pb * bg^2 + theta * bg
  M <- max(lw); M + log(sum(exp(lw - M)) * db)
}
Phi <- function(gam) {
  0.5 * lamg * (gam - mu0)^2 + 0.5 * P11re * gam^2 -
    sum(vapply(seq_len(J), function(j) A_j(Pb * gam, j), 0))
}

## ---- mean map M(gamma) ----
Mmap <- function(gam) {
  bbar <- vapply(seq_len(J), function(j) {
    lw <- lw_group(gam, j); w <- exp(lw - max(lw)); sum(w * bg) / sum(w)
  }, 0)
  (lamg * mu0 + Pb * sum(bbar)) / P11
}

gstar <- optimize(Phi, c(-6, 6), tol = 1e-10)$minimum
cat(sprintf("gamma*        = %.8f\n", gstar))
cat(sprintf("M(gamma*)     = %.8f   (fixed point residual %.2e)\n",
            Mmap(gstar), Mmap(gstar) - gstar))

## ---- q(gamma' | gamma) by 2-d quadrature over (beta_1, beta_2) ----
log_q <- function(gp, gam) {
  lw1 <- lw_group(gam, 1); lw2 <- lw_group(gam, 2)
  w1 <- exp(lw1 - max(lw1)); w2 <- exp(lw2 - max(lw2))
  w1 <- w1 / sum(w1); w2 <- w2 / sum(w2)
  mm <- outer(bg, bg, m_of)                        # m(beta_1, beta_2)
  dens <- dnorm(gp, mean = mm, sd = sqrt(1 / P11))
  log(sum(outer(w1, w2) * dens))
}
log_qQ <- function(gp) dnorm(gp, gstar, sqrt(1 / P11), log = TRUE)

Dfun <- function(gam) {
  optimize(function(gp) log_q(gp, gam) - log_qQ(gp), c(-12, 12), tol = 1e-9)$objective
}

## ---- checks ----
gs   <- seq(gstar - 1.6, gstar + 1.6, length.out = 13)
Dv   <- vapply(gs, Dfun, 0)
Dst  <- Dfun(gstar)
lhs  <- Dv - Dst
rhs  <- (vapply(gs, Phi, 0) - Phi(gstar)) - 0.5 * P11 * (gs - gstar)^2

cat("\n  gamma      D(g)-D(g*)     Phibar - quad        diff\n")
for (i in seq_along(gs))
  cat(sprintf("%7.3f  %13.8f  %13.8f  %11.2e\n", gs[i], lhs[i], rhs[i], lhs[i] - rhs[i]))

cat(sprintf("\nmax |identity residual|      = %.3e\n", max(abs(lhs - rhs))))

Psi <- Dst - Dv
cat(sprintf("min Psi (want >= 0)          = %.3e\n", min(Psi)))
cat(sprintf("max Psi - 0.5*P11*(g-g*)^2   = %.3e   (want <= 0)\n",
            max(Psi - 0.5 * P11 * (gs - gstar)^2)))
cat(sprintf("argmax epsilon               = %.6f  (gamma* = %.6f)\n",
            gs[which.max(Dv)], gstar))
d2 <- diff(diff(Psi)) / (gs[2] - gs[1])^2
cat(sprintf("min discrete Psi'' (want>=0) = %.3e\n", min(d2)))

## ---- coercivity: Psi''(gamma) = sum_j H_j' Pb V_j(gamma) Pb H_j, prior-free ----
## Cramer-Rao lower bound V_j >= (Pb + G_j)^{-1} with G_j = sup(-l_j'') = n_j/4 (logit)
S_lo <- sum(Pb^2 / (Pb + n_tr / 4))          # circumscribing metric
S_hi <- P11re                                 # = sum_j H' Pb H, Brascamp-Lieb ceiling
cat(sprintf("\nS (Cramer-Rao floor)         = %.6f\n", S_lo))
cat(sprintf("P11re (Brascamp-Lieb ceiling)= %.6f,  P11 = %.6f\n", S_hi, P11))
cat(sprintf("min Psi'' - S (want >= 0)    = %.3e\n", min(d2) - S_lo))
cat(sprintf("max Psi'' - P11re (want <=0) = %.3e\n", max(d2) - S_hi))

gg  <- seq(gstar - 5, gstar + 5, length.out = 21)
Pg  <- Dst - vapply(gg, Dfun, 0)
cat(sprintf("min Psi - 0.5*S*(g-g*)^2     = %.3e   (want >= 0: outer ellipsoid)\n",
            min(Pg - 0.5 * S_lo * (gg - gstar)^2)))
cat(sprintf("max Psi - 0.5*P11re*(g-g*)^2 = %.3e   (want <= 0: inner ellipsoid)\n",
            max(Pg - 0.5 * P11re * (gg - gstar)^2)))
cat(sprintf("Psi at |g-g*|=5              = %.4f  (coercive: grows)\n", Pg[length(Pg)]))
