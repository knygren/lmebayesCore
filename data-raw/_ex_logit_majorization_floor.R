# Worked logit example for inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md section 3A.
#
# Model: binomial random-intercept GLMM, scalar random effect, intercept-only
# hyper design, so q = 1 and every quantity below is a scalar.
#
#   y_j ~ Binomial(N_j, p(beta_j)),   beta_j ~ N(gamma, 1/lambda_b),
#   gamma ~ N(mu_0, 1/lambda_gamma)
#
# Group data precision, as in the note:   G_j(beta_j) = N_j p(beta_j)(1-p(beta_j)).
#
# Everything is computed by direct 1-D quadrature -- no sampler, no Laplace --
# so the "true" kappa is exact and the floor-based kappa^UB can be compared
# against it rather than against an approximation.
#
# Question asked: given a beta-escape budget delta_beta, how tight a box can we
# certify, what curvature floor does the box buy, and how much does capping the
# weights cost relative to the true conditional variances?
#
# Scratch example, not part of the package.

set.seed(11)

## ---- 1. simulate -----------------------------------------------------------

J        <- 20
N_j      <- c(rep(20, 6), rep(60, 7), rep(200, 7))       # trials per group
gamma_tr <- -1.0                                          # population log-odds
tau2     <- 0.5
lam_b    <- 1 / tau2                                      # RE precision, known
lam_g    <- 0.05                                          # weak population prior
mu_0     <- 0

beta_tr <- rnorm(J, gamma_tr, sqrt(tau2))
y_j     <- rbinom(J, N_j, plogis(beta_tr))

## ---- 2. grids --------------------------------------------------------------

bg <- seq(-9, 6, length.out = 4001)     # beta grid
gg <- seq(-4, 2, length.out = 601)      # gamma grid
db <- bg[2] - bg[1]; dg <- gg[2] - gg[1]

loglik <- sapply(1:J, function(j) y_j[j] * bg - N_j[j] * log1p(exp(bg)))  # 4001 x J

## ---- 3. marginal posterior of gamma, and gamma* ----------------------------

logZ <- sapply(1:J, function(j) {
  M <- outer(loglik[, j], rep(1, length(gg))) -
    0.5 * lam_b * outer(bg, gg, function(b, g) (b - g)^2)
  mx <- apply(M, 2, max)
  mx + log(colSums(exp(sweep(M, 2, mx))) * db)
})                                                        # 601 x J

lpost_g <- rowSums(logZ) - 0.5 * lam_g * (gg - mu_0)^2
w_g     <- exp(lpost_g - max(lpost_g)); w_g <- w_g / (sum(w_g) * dg)
gstar   <- gg[which.max(lpost_g)]

## refine gamma* by the mean map, which is what the note actually uses
bmean <- function(g) sapply(1:J, function(j) {
  lw <- loglik[, j] - 0.5 * lam_b * (bg - g)^2
  w  <- exp(lw - max(lw))
  sum(w * bg) / sum(w)
})
P11    <- lam_g + J * lam_b
P11_RE <- J * lam_b
for (it in 1:200) {
  gnew <- (lam_g * mu_0 + lam_b * sum(bmean(gstar))) / P11
  if (abs(gnew - gstar) < 1e-12) { gstar <- gnew; break }
  gstar <- gnew
}

## ---- 4. the TRUE kappa, from exact conditional variances at gamma* ---------

Vtrue <- sapply(1:J, function(j) {
  lw <- loglik[, j] - 0.5 * lam_b * (bg - gstar)^2
  w  <- exp(lw - max(lw)); w <- w / sum(w)
  m  <- sum(w * bg); sum(w * (bg - m)^2)
})
S_true   <- lam_b^2 * sum(Vtrue)
kap_true <- S_true / P11

## ---- 5. marginal posterior of each beta_j ----------------------------------

marg <- sapply(1:J, function(j) {
  M  <- outer(loglik[, j], rep(1, length(gg))) -
    0.5 * lam_b * outer(bg, gg, function(b, g) (b - g)^2)
  M  <- exp(sweep(M, 2, apply(M, 2, max)))
  M  <- sweep(M, 2, colSums(M) * db, "/")          # conditional densities
  as.vector(M %*% (w_g * dg))                      # integrate over gamma
})

## ---- 6. box for a given beta-budget, and the floor it buys -----------------

box_and_floor <- function(delta_beta) {
  per <- delta_beta / J
  out <- t(sapply(1:J, function(j) {
    cdf <- cumsum(marg[, j]) * db; cdf <- cdf / cdf[length(cdf)]
    lo  <- bg[which.max(cdf >= per / 2)]
    hi  <- bg[which.max(cdf >= 1 - per / 2)]
    inb <- bg >= lo & bg <= hi
    wmin <- min(plogis(bg[inb]) * (1 - plogis(bg[inb])))
    c(lo = lo, hi = hi, wmin = wmin, Gam = N_j[j] * wmin)
  }))
  data.frame(out)
}

certificate <- function(delta_beta) {
  bf   <- box_and_floor(delta_beta)
  omUB <- lam_b / (lam_b + bf$Gam)                 # capped pooling weights
  S_UB <- lam_b^2 * sum(1 / (lam_b + bf$Gam))      # = B^UB
  kUB  <- S_UB / P11
  PiLB <- P11 - S_UB
  nu   <- kUB / (1 - kUB)
  cprm <- (1 - kUB)^(-1 / 2)                       # q = 1
  list(bf = bf, om_max = max(omUB), kUB = kUB, PiLB = PiLB, nu = nu, cp = cprm)
}

## ---- 7. report -------------------------------------------------------------

cat("Simulated logit GLMM:  J =", J, " lambda_b =", lam_b,
    " lambda_gamma =", lam_g, "\n")
cat("gamma* =", format(gstar, digits = 4),
    "   true kappa =", format(kap_true, digits = 4),
    "   true rho = kappa/(1-kappa) =",
    format(kap_true / (1 - kap_true), digits = 4), "\n")
cat("group p at gamma*: ", paste(format(plogis(bmean(gstar)), digits = 2),
                                 collapse = " "), "\n\n")

cat("=== Trade-off: beta budget vs certified floor ===\n")
tab <- do.call(rbind, lapply(c(0.5, 0.2, 0.1, 0.05, 1e-2, 1e-3, 1e-4, 1e-6),
                             function(dB) {
  ce <- certificate(dB)
  data.frame(delta_beta = dB,
             half_width = mean((ce$bf$hi - ce$bf$lo) / 2),
             wmin_min = min(ce$bf$wmin),
             Gam_min  = min(ce$bf$Gam),
             om_max   = ce$om_max,
             kappa_UB = ce$kUB,
             rho      = ce$nu,
             PiLB_pos = ce$PiLB > 0,
             cprime   = ce$cp)
}))
print(tab, row.names = FALSE, digits = 3)

cat("\ntrue kappa =", format(kap_true, digits = 4),
    "-- the capped kappa_UB above is the price of the box.\n")

cat("\n=== Per-group floors at delta_beta = 1e-3 ===\n")
ce <- certificate(1e-3)
pg <- data.frame(j = 1:J, N = N_j, y = y_j,
                 p_hat = y_j / N_j,
                 lo = ce$bf$lo, hi = ce$bf$hi,
                 wmin = ce$bf$wmin, Gamma = ce$bf$Gam,
                 omega = lam_b / (lam_b + ce$bf$Gam))
print(pg[order(pg$omega, decreasing = TRUE), ], row.names = FALSE, digits = 3)

cat("\nomega is the pooling weight: the certificate is driven by the LARGEST,",
    "\ni.e. by the least informative group, not by the average.\n")

## ---- 8. does the floor degrade like exp(-C sqrt(log(1/delta)))? ------------
##
## Section 3A.4 predicts (a) box half-width linear in sqrt(log(1/delta_beta)),
## because the posterior tails are Gaussian, and (b) log(wmin) linear in the
## same, because the logit weight decays like exp(-|eta|).

cat("\n=== Section 3A.4 scaling check ===\n")
x <- sqrt(log(1 / tab$delta_beta))
f1 <- lm(tab$half_width ~ x)
f2 <- lm(log(tab$wmin_min) ~ x)
cat(sprintf("half_width ~ sqrt(log(1/delta)):  slope %.3f  R^2 %.4f\n",
            coef(f1)[2], summary(f1)$r.squared))
cat(sprintf("log(wmin)  ~ sqrt(log(1/delta)):  slope %.3f  R^2 %.4f\n",
            coef(f2)[2], summary(f2)$r.squared))
cat(sprintf("rho over six orders of magnitude in delta_beta: %.3f -> %.3f (x%.2f)\n",
            tab$rho[1], tab$rho[nrow(tab)], tab$rho[nrow(tab)] / tab$rho[1]))

## Honest comparison: a weak power law fits this range at least as well.  The
## sqrt(log) form is an ASYMPTOTIC prediction (Gaussian posterior tail composed
## with an exponentially decaying weight); over delta_beta in [0.5, 1e-6] the
## box half-width is only ~0.7-1.7, which is not deep enough into either tail
## for the two forms to separate.  Do not read the R^2 comparison as evidence
## for the asymptotic law.
f3 <- lm(log(tab$wmin_min) ~ log(1 / tab$delta_beta))
cat(sprintf("log(wmin)  ~ log(1/delta):        slope %.3f  R^2 %.4f\n",
            coef(f3)[2], summary(f3)$r.squared))
cat("both forms fit this range; they do not separate here.  The robust finding",
    "\nis the magnitude: rho roughly doubles over six orders of magnitude in",
    "\ndelta_beta, whichever functional form is used to describe it.\n")
