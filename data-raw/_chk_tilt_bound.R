# Numerical checks for inst/CHAPTER_C05_IMPLEMENTATION.md section 4B.
#
# Check A (4B.1/4B.3 algebra). lambda_star = 1/k0 - q/(2d) minimizes
#   exp(-(lam-1) d) (1 - lam k0)^(-q/2)   over lam in (1, 1/k0),
# the minimum equals the closed form (2 e d / (q k0))^(q/2) exp(-d (1-k0)/k0),
# and the side condition lambda_star > 1 is exactly d > q k0 / (2 (1-k0)).
#
# Check B (4B.2 MGF). Under closure, Psi = (1/2) sum k_i Z_i^2 with Z ~ N(0,I)
# under Q, so M_Q(lam) = prod (1 - lam k_i)^(-1/2) for lam < 1/max(k_i).
#
# Check C (4B.0 tilting identity). Under closure,
#   E_Q[e^Psi 1{Psi > d}] / E_Q[e^Psi]  =  Pr(sum mu_i Z_i^2 > 2d),
# mu_i = k_i / (1 - k_i) -- the exact escape probability of section 4A.3.
#
# Check D (4B.3 domination). The closed form dominates the exact escape
# probability whenever the side condition holds, for flat and non-flat spectra.
#
# Scratch check, not part of the package.

set.seed(20260816)

cat("=== Check A: lambda_star algebra and closed form ===\n")

grid <- expand.grid(
  q  = c(1, 2, 3, 5, 8),
  k0 = c(0.2, 0.5, 0.8, 0.95),
  d  = c(2, 5, 10, 30, 100)
)

outA <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  q <- grid$q[i]; k0 <- grid$k0[i]; d <- grid$d[i]

  logbound <- function(lam) -(lam - 1) * d - (q / 2) * log(1 - lam * k0)

  lam_star <- 1 / k0 - q / (2 * d)
  side_ok  <- lam_star > 1
  side_alg <- d > q * k0 / (2 * (1 - k0))

  # logbound is convex on (1, 1/k0) and blows up at the upper end, so the
  # infimum is interior exactly when the side condition holds.
  num <- optimize(logbound, c(1 + 1e-10, 1 / k0 - 1e-12))

  closed <- (2 * exp(1) * d / (q * k0))^(q / 2) * exp(-d * (1 - k0) / k0)
  exact  <- pchisq(2 * d * (1 - k0) / k0, df = q, lower.tail = FALSE)

  data.frame(
    q = q, k0 = k0, d = d,
    side_ok = side_ok,
    side_matches_algebra = identical(side_ok, side_alg),
    lam_gap = abs(lam_star - num$minimum),
    rel_gap = abs(closed - exp(num$objective)) / exp(num$objective),
    closed = closed, exact = exact,
    dominates = closed >= exact
  )
}))

keepA <- outA[outA$side_ok, ]

cat("cases total / meeting side condition:",
    nrow(outA), "/", nrow(keepA), "\n")
cat("side condition matches d > q k0 / (2 (1-k0)):",
    all(outA$side_matches_algebra), "\n")
cat("max |lam_star - numeric argmin|:      ",
    format(max(keepA$lam_gap), digits = 3), "\n")
cat("max rel gap closed vs numeric inf:    ",
    format(max(keepA$rel_gap), digits = 3), "\n")

cat("\n=== Check B: closure MGF M_Q(lam) = prod (1 - lam k_i)^(-1/2) ===\n")

n  <- 4e6
kB <- c(0.7, 0.4, 0.15)
Z  <- matrix(rnorm(n * length(kB)), ncol = length(kB))
Psi <- 0.5 * as.vector(Z^2 %*% kB)

for (lam in c(0.5, 1, 1.2)) {
  mc     <- mean(exp(lam * Psi))
  closed <- prod((1 - lam * kB)^(-1 / 2))
  cat(sprintf("lam = %4.2f   MC = %10.6f   closed = %10.6f   rel err = %.2e\n",
              lam, mc, closed, abs(mc - closed) / closed))
}

cat("\n=== Check C: tilting identity reproduces the exact escape probability ===\n")

muB   <- kB / (1 - kB)
denom <- mean(exp(Psi))
for (d in c(1, 2, 4, 6)) {
  tilt  <- mean(exp(Psi) * (Psi > d)) / denom
  exact <- mean(0.5 * as.vector(matrix(rnorm(n * length(muB)),
                                       ncol = length(muB))^2 %*% muB) > d)
  cat(sprintf("d = %d   tilt identity = %.6f   direct mu-draw = %.6f   rel err = %.2e\n",
              d, tilt, exact, abs(tilt - exact) / exact))
}

cat("\n=== Check D: domination, flat and non-flat spectra ===\n")

cat("flat spectrum, all side-condition cases dominate:",
    all(keepA$dominates), " (violations:", sum(!keepA$dominates), ")\n")

nonflat <- do.call(rbind, lapply(1:400, function(i) {
  q  <- sample(2:8, 1)
  k0 <- runif(1, 0.15, 0.97)
  ks <- sort(c(k0, runif(q - 1, 0.01, k0)), decreasing = TRUE)
  d  <- runif(1, q * k0 / (2 * (1 - k0)) * 1.01, 200)
  closed <- (2 * exp(1) * d / (q * k0))^(q / 2) * exp(-d * (1 - k0) / k0)
  mu     <- ks / (1 - ks)
  ex     <- mean(0.5 * as.vector(matrix(rnorm(2e5 * q), ncol = q)^2 %*% mu) > d)
  data.frame(q = q, k0 = k0, d = d, closed = closed, exact = ex,
             dominates = closed >= ex)
}))

cat("non-flat spectra, cases dominating:",
    sum(nonflat$dominates), "/", nrow(nonflat), "\n")

cat("\nconservatism ratio closed/exact, flat spectrum:\n")
keepA$ratio <- keepA$closed / keepA$exact
print(summary(keepA$ratio))

cat("\n=== Check E: what the bound is worth at Prior_Setup_GLMM's pop.pwt ===\n")
#
# Under the Prior_Setup convention Sigma_fixef = ((1-w)/w) vcov(fit_ref), so
# Lambda_gamma = (w/(1-w)) vcov(fit_ref)^{-1}.  If vcov(fit_ref)^{-1} aligns with
# P11_RE then P11 = P11_RE/(1-w) and 1 - k0 = Lambda_gamma / P11 = w exactly, so
# k0 = 1 - pop.pwt and the decay exponent is w/(1-w) = pop.nprior / J.
#
# Certificate arithmetic: eps = eps(gamma*) e^{-d} Q(Ctilde_d), with the closure
# floor eps(gamma*) >= (1 + k0)^{-q/2} and Q(Ctilde_d) >= pchisq(2d/k0, q).

q <- 3; delta <- 0.01; tol <- 0.05

d_for_delta <- function(k0, q, delta) {
  bnd <- function(d) (2 * exp(1) * d / (q * k0))^(q / 2) * exp(-d * (1 - k0) / k0)
  lo <- q * k0 / (2 * (1 - k0)) * 1.0001
  hi <- lo
  while (bnd(hi) > delta && hi < 1e12) hi <- hi * 2
  if (hi >= 1e12) return(NA_real_)
  uniroot(function(d) bnd(d) - delta, c(lo, hi), tol = 1e-8)$root
}

tabE <- do.call(rbind, lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9),
                              function(w) {
  k0 <- 1 - w
  d  <- d_for_delta(k0, q, delta)
  eps_star <- (1 + k0)^(-q / 2)
  eps <- eps_star * exp(-d) * pchisq(2 * d / k0, df = q)
  n   <- if (eps > 0) log(tol - delta) / log1p(-eps) else Inf
  data.frame(pop.pwt = w, nprior_over_J = w / (1 - w), k0 = k0,
             d = d, eps = eps, n_sweeps = n)
}))

print(tabE, row.names = FALSE, digits = 4)

cat("\nInterpretation: exponent (1-k0)/k0 = pop.pwt/(1-pop.pwt) = pop.nprior/J.\n")
cat("Package defaults are pop.pwt_default_low = 0.01, high = 0.05.\n")
