# Verification of the closed forms in inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md
# section 3A.4: the Hessian sandwich (S'), the outer-ellipsoid floor (F), and the
# level-set mass bound (D).
#
# Same model as data-raw/_ex_logit_majorization_floor.R (same seed and design):
# binomial random-intercept logit, scalar RE, intercept-only hyper design, one
# binomial observation of N_j trials per group, so z_j = 1 and q = 1.
#
# The marginal after integrating gamma out of the PRIOR (LOGIT_MARGINAL_INTEGRATE_GAMMA.md
# section 5, with mu_0 = 0, x_j = 1) is
#
#   pi~(beta | y)  =  exp(sum_j l_j(beta_j)) * exp(-0.5 * beta' Lam_b beta),
#   Lam_b = lam_b I - lam_b^2/(lam_g + J lam_b) 11',
#
# and Xi(beta) = log pi~(beta_dag) - log pi~(beta), B_r = {Xi <= r}.
#
# Exact iid draws from pi(beta | y) are available WITHOUT a sampler: draw gamma
# from its 1-D grid posterior, then the beta_j independently from their 1-D
# conditionals.  That is what check (D) is compared against.
#
# Scratch check, not part of the package.

set.seed(11)

## ---- model, identical to the worked example --------------------------------

J        <- 20
N_j      <- c(rep(20, 6), rep(60, 7), rep(200, 7))
gamma_tr <- -1.0
tau2     <- 0.5
lam_b    <- 1 / tau2
lam_g    <- 0.05
mu_0     <- 0

beta_tr <- rnorm(J, gamma_tr, sqrt(tau2))
y_j     <- rbinom(J, N_j, plogis(beta_tr))

P11_pri <- lam_g + J * lam_b                       # prior P_11 (x_j = 1)
Lam_b   <- diag(lam_b, J) - (lam_b^2 / P11_pri) * matrix(1, J, J)
Lam_UB  <- Lam_b + diag(N_j / 4)                   # global logit ceiling w <= n/4

stopifnot(min(eigen(Lam_b, symmetric = TRUE)$values) > 0)

## ---- Xi, its gradient and Hessian ------------------------------------------

negll  <- function(b) -sum(y_j * b - N_j * log1p(exp(b)))
gr_nll <- function(b) -(y_j - N_j * plogis(b))
He_nll <- function(b) diag(N_j * plogis(b) * (1 - plogis(b)))

f      <- function(b) negll(b) + 0.5 * drop(crossprod(b, Lam_b %*% b))
gr_f   <- function(b) gr_nll(b) + drop(Lam_b %*% b)
He_f   <- function(b) He_nll(b) + Lam_b

newton <- function(b, tgt = function(x) gr_f(x)) {
  for (i in 1:200) {
    g <- tgt(b); s <- solve(He_f(b), g); b <- b - s
    if (max(abs(s)) < 1e-12) break
  }
  b
}

b_dag <- newton(rep(0, J))
f_dag <- f(b_dag)
Xi    <- function(b) f(b) - f_dag          # normalizer-free by construction

cat("mode found, max|grad| =", format(max(abs(gr_f(b_dag))), digits = 3), "\n\n")

## ---- Check 1: the sandwich (S') --------------------------------------------

cat("=== Check 1: sandwich (S'), 0.5||.||^2_Lam_b <= Xi <= 0.5||.||^2_Lam_UB ===\n")
qform <- function(b, A) drop(crossprod(b - b_dag, A %*% (b - b_dag)))
viol_lo <- viol_hi <- numeric(0)
for (k in 1:4000) {
  sc <- 10^runif(1, -1.5, 0.8)                      # a wide range of radii
  b  <- b_dag + sc * rnorm(J)
  viol_lo <- c(viol_lo, Xi(b) - 0.5 * qform(b, Lam_b))
  viol_hi <- c(viol_hi, 0.5 * qform(b, Lam_UB) - Xi(b))
}
cat(sprintf("lower gap: min %.3e  (must be >= 0)   %s\n", min(viol_lo),
            if (min(viol_lo) >= -1e-9) "OK" else "FAIL"))
cat(sprintf("upper gap: min %.3e  (must be >= 0)   %s\n", min(viol_hi),
            if (min(viol_hi) >= -1e-9) "OK" else "FAIL"))

## ---- Check 2: (F), outer ellipsoid vs the true support function ------------
##
## True max of beta_j over {Xi <= r}: minimise Xi(beta) - t*beta_j for t > 0,
## then choose t so the solution lands on the level surface.

cat("\n=== Check 2: (F) outer-ellipsoid bound on max |eta_j| over B_r ===\n")
Lbi <- solve(Lam_b)

sup_coord <- function(r, j, sgn) {
  e <- numeric(J); e[j] <- sgn
  g_t <- function(t) {
    b <- newton(b_dag, function(x) gr_f(x) - t * e)
    list(b = b, xi = Xi(b))
  }
  hi <- 1
  while (g_t(hi)$xi < r && hi < 1e8) hi <- hi * 2
  t <- uniroot(function(t) g_t(t)$xi - r, c(1e-10, hi), tol = 1e-10)$root
  drop(g_t(t)$b[j])
}

res2 <- do.call(rbind, lapply(c(2, 8, 20), function(r) {
  do.call(rbind, lapply(c(3, 7, 18), function(j) {
    tru <- max(abs(c(sup_coord(r, j, +1), sup_coord(r, j, -1))))
    bnd <- abs(b_dag[j]) + sqrt(2 * r * Lbi[j, j])       # formula (F)
    data.frame(r = r, j = j, N = N_j[j], true_max = tru, bound_F = bnd,
               ok = bnd >= tru - 1e-8, ratio = bnd / tru)
  }))
}))
print(res2, row.names = FALSE, digits = 4)
cat("all (F) bounds valid:", all(res2$ok),
    "   median conservatism factor:", format(median(res2$ratio), digits = 3), "\n")

## ---- Check 3: (D), level-set mass bound vs exact iid draws -----------------

cat("\n=== Check 3: (D) mass bound vs exact draws from pi(beta|y) ===\n")

## eigen identity first
th   <- eigen(solve(Lam_b) %*% Lam_UB, only.values = TRUE)$values
th   <- sort(Re(th))
cat(sprintf("prod(theta) = %.6e   det ratio = %.6e   rel err %.2e\n",
            prod(th), det(Lam_UB) / det(Lam_b),
            abs(prod(th) - det(Lam_UB) / det(Lam_b)) / (det(Lam_UB) / det(Lam_b))))

## exact iid draws: gamma ~ pi(gamma|y) on a grid, then beta_j | gamma, y
bg <- seq(-9, 6, length.out = 3001); db <- bg[2] - bg[1]
gg <- seq(-4, 2, length.out = 401)
loglik <- sapply(1:J, function(j) y_j[j] * bg - N_j[j] * log1p(exp(bg)))

logZ <- sapply(1:J, function(j) {
  M  <- outer(loglik[, j], rep(1, length(gg))) -
    0.5 * lam_b * outer(bg, gg, function(b, g) (b - g)^2)
  mx <- apply(M, 2, max)
  mx + log(colSums(exp(sweep(M, 2, mx))) * db)
})
lpg <- rowSums(logZ) - 0.5 * lam_g * (gg - mu_0)^2
wg  <- exp(lpg - max(lpg)); wg <- wg / sum(wg)

nmc  <- 2e5L
cnt  <- drop(rmultinom(1, nmc, wg))
draw <- matrix(0, nmc, J); pos <- 0L
for (gi in which(cnt > 0)) {
  n_g <- cnt[gi]
  for (j in 1:J) {
    lw <- loglik[, j] - 0.5 * lam_b * (bg - gg[gi])^2
    w  <- exp(lw - max(lw)); cdf <- cumsum(w) / sum(w)
    draw[(pos + 1L):(pos + n_g), j] <-
      bg[findInterval(runif(n_g), cdf) + 1L]
  }
  pos <- pos + n_g
}
draw <- draw[sample(nmc), , drop = FALSE]

## sanity: the draws must reproduce the grid marginal of one coordinate
mg3 <- {
  M <- outer(loglik[, 3], rep(1, length(gg))) -
    0.5 * lam_b * outer(bg, gg, function(b, g) (b - g)^2)
  M <- exp(sweep(M, 2, apply(M, 2, max)))
  M <- sweep(M, 2, colSums(M) * db, "/")
  as.vector(M %*% wg)
}
cat(sprintf("sanity: E[beta_3] grid %.4f vs draws %.4f;  sd %.4f vs %.4f\n",
            sum(mg3 * bg) * db, mean(draw[, 3]),
            sqrt(sum(mg3 * (bg - sum(mg3 * bg) * db)^2) * db), sd(draw[, 3])))

Xi_draw <- apply(draw, 1, Xi)

wchisq_upper <- function(t, w, nmc = 4e5L) {
  Z <- matrix(rnorm(nmc * length(w)), nmc, length(w))
  mean(drop(Z^2 %*% w) > t)
}
set.seed(99)
res3 <- do.call(rbind, lapply(c(10, 15, 20, 25, 30, 40), function(r) {
  exact <- mean(Xi_draw > r)
  bnd   <- sqrt(prod(th)) * wchisq_upper(2 * r, th)
  data.frame(r = r, exact_mass = exact, bound_D = min(1, bnd),
             ok = bnd >= exact - 1e-12, ratio = bnd / max(exact, 1e-12))
}))
print(res3, row.names = FALSE, digits = 4)
cat("all (D) bounds valid:", all(res3$ok), "\n")
cat("constant sqrt(prod(theta)) =", format(sqrt(prod(th)), digits = 4), "\n")

## ---- Check 4: does the bootstrap rescue (F)? -------------------------------
##
## Iterate  Lam_LB <- Lam_b + diag(Gamma_j),  Gamma_j = N_j * w(eta_max_j),
## eta_max_j = |b_dag_j| + sqrt(2 r (Lam_LB)^{-1}_{jj}).
##
## Validity: if grad^2 Xi >= Lam_LB on E = {||.||^2_Lam_LB <= 2r}, then for any
## beta outside E the segment from b_dag to beta leaves E at some beta_bar with
## Xi(beta_bar) >= r, and Xi increases along rays from the mode (convexity plus
## grad Xi(b_dag) = 0), so Xi(beta) >= r.  Hence B_r is contained in E, which is
## what (F) needs.  The iteration is monotone, so it converges.

cat("\n=== Check 4: bootstrapped outer ellipsoid vs plain (F) ===\n")
## The map T(L) = Lam_b + diag(N_j w(|b_dag| + sqrt(2 r (L^{-1})_jj))) is MONOTONE
## increasing in L: larger L shrinks the ellipsoid, raises w, raises T(L).  By
## Knaster-Tarski its fixed points form a lattice, and iterating from below
## (L = Lam_b) reaches the LEAST fixed point -- which is the trivial one whenever
## the starting ellipsoid is so wide that w ~ 0.  Iterating from above reaches
## the GREATEST fixed point.  Only the fixed point itself has to satisfy the
## validity condition, not the iterates, so starting high is legitimate.
boot_eta <- function(r, from = c("below", "above"), iters = 500) {
  from <- match.arg(from)
  L <- if (from == "below") Lam_b else Lam_b + diag(N_j / 4)
  for (k in 1:iters) {
    em <- abs(b_dag) + sqrt(2 * r * diag(solve(L)))
    Ln <- Lam_b + diag(N_j * plogis(em) * (1 - plogis(em)))
    if (max(abs(Ln - L)) < 1e-12) { L <- Ln; break }
    L <- Ln
  }
  list(eta = abs(b_dag) + sqrt(2 * r * diag(solve(L))), L = L)
}

res4 <- do.call(rbind, lapply(c(2, 8, 20, 40), function(r) {
  lo <- boot_eta(r, "below"); hi <- boot_eta(r, "above")
  do.call(rbind, lapply(c(3, 7, 18), function(j) {
    tru <- max(abs(c(sup_coord(r, j, +1), sup_coord(r, j, -1))))
    data.frame(r = r, j = j, N = N_j[j], true_max = tru,
               plain_F = abs(b_dag[j]) + sqrt(2 * r * Lbi[j, j]),
               from_below = lo$eta[j], from_above = hi$eta[j],
               ok = hi$eta[j] >= tru - 1e-8,
               ratio = hi$eta[j] / tru)
  }))
}))
print(res4, row.names = FALSE, digits = 4)
cat("greatest fixed point still valid:", all(res4$ok),
    "  median conservatism: plain",
    format(median(res2$ratio), digits = 3),
    "-> from above", format(median(res4$ratio), digits = 3), "\n")

cat("\nImplied weight floors at r = 8 (true vs plain vs greatest fixed point):\n")
bb8 <- boot_eta(8, "above")
tru8 <- sapply(1:J, function(j)
  max(abs(c(sup_coord(8, j, +1), sup_coord(8, j, -1)))))
wof <- function(e) plogis(e) * (1 - plogis(e))
cmp <- data.frame(j = 1:J, N = N_j,
                  w_true = wof(tru8),
                  w_plain = wof(abs(b_dag) + sqrt(2 * 8 * diag(Lbi))),
                  w_boot = wof(bb8$eta))
cmp$om_true  <- lam_b / (lam_b + N_j * cmp$w_true)
cmp$om_plain <- lam_b / (lam_b + N_j * cmp$w_plain)
cmp$om_boot  <- lam_b / (lam_b + N_j * cmp$w_boot)
print(cmp, row.names = FALSE, digits = 3)
cat(sprintf("\nomega_max: true %.4f   plain(F) %.6f   bootstrapped %.4f\n",
            max(cmp$om_true), max(cmp$om_plain), max(cmp$om_boot)))

## ---- Check 5: the phase transition in r ------------------------------------
##
## The greatest fixed point is nontrivial only while the level set is small
## enough for the weights to hold the ellipsoid in.  Past a critical r* the only
## fixed point is Lam_b itself and the certificate goes vacuous -- the Theorem 2
## trade-off appearing as a sharp threshold rather than a smooth decay.

cat("\n=== Check 5: omega_max against r (greatest fixed point) ===\n")
res5 <- do.call(rbind, lapply(c(2, 5, 8, 12, 16, 20, 24, 26, 28, 30, 35, 40),
                              function(r) {
  e <- boot_eta(r, "above")$eta
  om <- lam_b / (lam_b + N_j * wof(e))
  data.frame(r = r, delta_beta = mean(Xi_draw > r), eta_max = max(e),
             wmin = min(wof(e)), om_max = max(om), trivial = max(e) > 20)
}))
print(res5, row.names = FALSE, digits = 4)
r_star <- max(res5$r[!res5$trivial])
cat(sprintf("\nlargest r with a nontrivial floor: between %.0f and %.0f\n",
            r_star, min(res5$r[res5$trivial])))
cat("Exact-draw mass outside B_r at that level:",
    format(mean(Xi_draw > r_star), digits = 3), "\n")
cat("So the certificate is capped at a beta-budget of roughly this size;",
    "\ntightening delta_beta below it cannot be done by this route.\n")

## ---- Check 6: level set vs box, head to head at matched budget -------------
##
## Decisive comparison, using the SAME exact draws for both shapes and the TRUE
## support function for the level set (no ellipsoid relaxation), so the contest
## is purely about shape.
##
##   box       B = prod_j [lo_j, hi_j], equal-tailed, tail mass delta/J per group
##   level set B_r with r the (1-delta) quantile of Xi
##
## Both carry mass >= 1 - delta.  Compare omega_max = max_j lam_b/(lam_b+Gamma_j),
## which is the only feature of the restriction the certificate depends on.

cat("\n=== Check 6: level set vs box at matched delta_beta (both exact) ===\n")
res6 <- do.call(rbind, lapply(c(0.5, 0.2, 0.05, 0.01), function(dB) {
  ## box
  qs <- sapply(1:J, function(j)
    quantile(draw[, j], c(dB / J / 2, 1 - dB / J / 2), names = FALSE))
  w_box <- pmin(wof(abs(qs[1, ])), wof(abs(qs[2, ])))
  om_box <- max(lam_b / (lam_b + N_j * w_box))
  ## level set, true support function
  r  <- unname(quantile(Xi_draw, 1 - dB))
  em <- sapply(1:J, function(j)
    max(abs(c(sup_coord(r, j, +1), sup_coord(r, j, -1)))))
  om_ls <- max(lam_b / (lam_b + N_j * wof(em)))
  data.frame(delta_beta = dB, r = r,
             box_maxeta = max(abs(qs)), ls_maxeta = max(em),
             om_box = om_box, om_levelset = om_ls,
             winner = if (om_box < om_ls) "box" else "level set")
}))
print(res6, row.names = FALSE, digits = 4)
cat("\nThe floor depends on the WORST-CASE coordinate excursion.  A superlevel",
    "\nset of a sum lets one coordinate run far while the others sit near the",
    "\nmode; a product set forbids exactly that.  Hence the box wins on this",
    "\nfunctional even though the level set is the better-shaped set overall.\n")

## ---- Check 7: the RAY mass bound, replacing (D) -----------------------------
##
## Polar coordinates at the mode.  Along the ray beta_dag + t u write
## phi_u(t) = Xi(beta_dag + t u); it is convex with phi_u(0) = 0, so phi_u(t)/t
## is nondecreasing.  With t_r(u) the crossing of the level surface,
##
##   t >= t_r :  phi_u(t) >= r t / t_r      (used past the boundary)
##   t <= t_r :  phi_u(t) <= r t / t_r      (used inside)
##
## Both radial integrals then contribute the SAME factor (t_r(u)/r)^n, which
## cancels in the ratio, leaving a bound that depends only on n and r:
##
##   pi(B_r^c | y)  <=  Gamma(n, r) / gamma(n, r).
##
## This is the fix for the vacuous constant in (D): past the boundary the
## curvature is DECLINING along every ray, so the density decays exponentially
## in t, not like a Gaussian -- which is exactly why a Gaussian majorant with
## prior precision needed a det-ratio constant exponential in n.

cat("\n=== Check 7: ray mass bound Gamma(n,r)/gamma(n,r) vs exact ===\n")
n_dim <- J
ray_bound <- function(r, n = n_dim)
  pgamma(r, n, lower.tail = FALSE) / pgamma(r, n, lower.tail = TRUE)

res7 <- do.call(rbind, lapply(c(10, 15, 20, 25, 30, 35, 40, 50), function(r) {
  ex <- mean(Xi_draw > r)
  data.frame(r = r, exact = ex, ray = min(1, ray_bound(r)),
             old_D = min(1, sqrt(prod(th)) * wchisq_upper(2 * r, th)),
             ok = ray_bound(r) >= ex - 1e-9,
             ratio = ray_bound(r) / max(ex, 1e-9))
}))
print(res7, row.names = FALSE, digits = 4)
cat("ray bound valid everywhere:", all(res7$ok), "\n")

## direct check of the radial convexity inequalities the proof rests on
bad <- 0
for (k in 1:300) {
  u  <- rnorm(J); u <- u / sqrt(sum(u^2))
  ph <- function(t) Xi(b_dag + t * u)
  r0 <- 12
  hi <- 1; while (ph(hi) < r0) hi <- hi * 2
  tr <- uniroot(function(t) ph(t) - r0, c(1e-8, hi), tol = 1e-12)$root
  ins <- runif(5, 0, tr); out <- tr * (1 + rexp(5))
  if (any(sapply(ins, ph) > r0 * ins / tr + 1e-8)) bad <- bad + 1
  if (any(sapply(out, ph) < r0 * out / tr - 1e-8)) bad <- bad + 1
}
cat("radial chord inequalities violated in", bad, "of 300 random directions\n")

## ---- Check 8: is the Check 6 comparison even admissible? -------------------
##
## Check 6 sized the box from EXACT per-group marginals pi(beta_j | y).  For
## J > 1 those have no closed form (LOGIT_MARGINAL_INTEGRATE_GAMMA.md section 8),
## and here they were obtained only by grid-integrating gamma out -- i.e. by the
## very computation the certificate is supposed to avoid.  So the "box" of
## Check 6 is an ORACLE, not a constructible set.
##
## The constructible box is the projection box of the level set:
##   proj_j = [ min_{B_r} beta_j , max_{B_r} beta_j ]   (support function)
##   PB     = prod_j proj_j  >=  B_r,   so  pi(PB^c) <= pi(B_r^c) <= (P2).
## Its omega_max is IDENTICAL to the level set's, since omega depends on the set
## only through its coordinate projections.  So level set and constructible box
## are the same certificate, and the oracle box is not a rival construction.

cat("\n=== Check 8: oracle box vs constructible (projection) box ===\n")
res_adm <- do.call(rbind, lapply(c(0.5, 0.2, 0.05), function(dB) {
  r <- uniroot(function(r) ray_bound(r) - dB, c(1, 200), tol = 1e-10)$root
  em <- sapply(1:J, function(j)
    max(abs(c(sup_coord(r, j, +1), sup_coord(r, j, -1)))))
  om_ls <- max(lam_b / (lam_b + N_j * wof(em)))
  inbox <- rowSums(sweep(abs(draw), 2, em, ">")) == 0
  qs <- sapply(1:J, function(j)
    quantile(draw[, j], c(dB / J / 2, 1 - dB / J / 2), names = FALSE))
  om_or <- max(lam_b / (lam_b + N_j * pmin(wof(abs(qs[1, ])), wof(abs(qs[2, ])))))
  data.frame(certified_delta = dB, r = r,
             om_levelset_certified = om_ls,
             om_oracle_box = om_or,
             exact_mass_levelset = mean(Xi_draw > r),
             exact_mass_projbox = mean(!inbox))
}))
print(res_adm, row.names = FALSE, digits = 4)
cat("\nThe oracle box needs per-group marginals, which are exactly what is not",
    "\navailable without sampling.  The projection box is constructible and has",
    "\nthe same omega_max as the level set, with a certified budget from (P2).\n")

cat("\nLevel r needed for a given certified budget, and the floor there:\n")
res8 <- do.call(rbind, lapply(c(0.5, 0.1, 0.05, 0.01), function(dB) {
  r <- uniroot(function(r) ray_bound(r) - dB, c(1, 200), tol = 1e-10)$root
  e <- boot_eta(r, "above")$eta
  data.frame(certified_delta = dB, r = r, exact_mass = mean(Xi_draw > r),
             om_max_boot = max(lam_b / (lam_b + N_j * wof(e))))
}))
print(res8, row.names = FALSE, digits = 4)
