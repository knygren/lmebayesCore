# The map beta_j -> G_j(beta_j), and the matrix floor it supports, for p_re = 2.
# Companion to section 3A.4 of inst/GAUSSIAN_MAJORIZATION_ESCAPE_BOUND.md.
#
# THE MAP (binomial logit, group j, observations i = 1..n_j):
#
#   eta_{ij}(beta_j) = o_{ij} + z_{ij}' beta_j                  (linear, R^p_re -> R)
#   p_{ij}           = expit(eta_{ij})
#   w_{ij}(beta_j)   = n_{ij} p_{ij} (1 - p_{ij})               (scalar weight)
#   G_j(beta_j)      = sum_i w_{ij}(beta_j) z_{ij} z_{ij}'      (p_re x p_re)
#                    = Z_j' W_j(beta_j) Z_j,  W_j = diag(w_{ij})
#
# Two structural facts drive everything:
#
#   (i)  w_{ij} depends on beta_j ONLY through the single linear functional
#        eta_{ij}.  So its minimum over a convex set is a one-dimensional
#        question, whatever p_re is.
#   (ii) log w(eta) = const + eta - 2 log(1 + e^eta) has second derivative
#        -2 p(1-p) < 0, so w is STRICTLY LOG-CONCAVE in eta, peaks at eta = 0,
#        and is strictly decreasing in |eta|.  Hence on an interval [a,b],
#              min w = min( w(a), w(b) )
#        with no interior case to consider.
#
# Consequently the weights are found from 2 n_j support-function evaluations per
# group REGARDLESS of p_re; only the assembly of G_j changes with dimension.
#
# Floors, both valid because w_i >= wbar_i pointwise and z z' >= 0:
#   (1) scalar   Gamma^1 = min_i wbar_{ij} * Z_j'Z_j
#   (2) per-obs  Gamma^2 = sum_i wbar_{ij} z_{ij} z_{ij}'        <- recommended
# Neither is attained at any single beta_j: different observations reach their
# minima at different points.  That is fine -- validity is all that is needed,
# and it is exactly the same "no joint minimiser required" argument that lets
# the floor be assembled group by group.
#
# (A third, tighter option -- the largest Gamma with v'Gamma v <= min_{B_j}
# v'G_j v for every unit v -- is a small SDP in 2-D.  It is omitted here: the
# characterisation is correct, but enforcing it on a finite grid of directions
# and evaluating the inner minimum on a finite sample both break feasibility,
# and a discretised version tested earlier produced a Gamma that failed the PSD
# check.  It is a target, not a certificate.)
#
# Scratch check, not part of the package.

set.seed(7)

J <- 6; pre <- 2; q <- 2; nobs <- 8
n_tr  <- 5
lam_g <- 0.05
Pb    <- matrix(c(2.0, 0.3, 0.3, 3.0), 2, 2)
mu0   <- c(-0.6, 0.4)

Zl <- lapply(1:J, function(j) cbind(1, rnorm(nobs, 0, 1)))
b_true <- t(sapply(1:J, function(j) mu0 + drop(backsolve(chol(Pb), rnorm(2)))))
yl <- lapply(1:J, function(j)
  rbinom(nobs, n_tr, plogis(drop(Zl[[j]] %*% b_true[j, ]))))

n_dim <- J * pre
Wj    <- diag(pre)
idx   <- function(j) ((j - 1) * pre + 1):(j * pre)

P11p  <- lam_g * diag(q) + J * t(Wj) %*% Pb %*% Wj
Lam_b <- matrix(0, n_dim, n_dim)
for (j in 1:J) for (k in 1:J)
  Lam_b[idx(j), idx(k)] <- (if (j == k) Pb else matrix(0, pre, pre)) -
    Pb %*% Wj %*% solve(P11p, t(Wj) %*% Pb)
Lam_b <- (Lam_b + t(Lam_b)) / 2
mu_b  <- rep(drop(Wj %*% mu0), J)
stopifnot(min(eigen(Lam_b, symmetric = TRUE)$values) > 1e-10)

## ---- the map, coded exactly as written above -------------------------------

eta_j <- function(j, bj) drop(Zl[[j]] %*% bj)
w_j   <- function(j, bj) { p <- plogis(eta_j(j, bj)); n_tr * p * (1 - p) }
G_j   <- function(j, bj) {                       # sum form
  Z <- Zl[[j]]; w <- w_j(j, bj)
  Reduce(`+`, lapply(1:nobs, function(i) w[i] * tcrossprod(Z[i, ])))
}
G_j_mat <- function(j, bj) t(Zl[[j]]) %*% (w_j(j, bj) * Zl[[j]])   # Z'WZ form

cat("=== The map: sum_i w_i z_i z_i'  ==  Z' W Z ? ===\n")
d <- max(sapply(1:20, function(k) {
  j <- sample(J, 1); bj <- rnorm(pre, 0, 1.5)
  max(abs(G_j(j, bj) - G_j_mat(j, bj)))
}))
cat("max discrepancy over 20 random (j, beta_j):", format(d, digits = 3), "\n")

bj0 <- b_true[1, ]
cat("\nWorked instance, group 1 at beta_1 =", format(bj0, digits = 3), ":\n")
print(data.frame(i = 1:nobs, z1 = Zl[[1]][, 1], z2 = round(Zl[[1]][, 2], 3),
                 eta = round(eta_j(1, bj0), 3),
                 p = round(plogis(eta_j(1, bj0)), 3),
                 w = round(w_j(1, bj0), 4)), row.names = FALSE)
cat("G_1(beta_1) =\n"); print(round(G_j(1, bj0), 4))

## ---- fact (ii): w is log-concave in eta, so minima sit at interval endpoints -

cat("\n=== w(eta) = n p(1-p): log-concavity and endpoint minima ===\n")
ee <- seq(-8, 8, length.out = 20001)
lw <- log(n_tr * plogis(ee) * (1 - plogis(ee)))
d2 <- diff(lw, differences = 2)
cat("max second difference of log w (must be < 0):", format(max(d2), digits = 3), "\n")
bad <- 0
for (k in 1:2000) {
  ab <- sort(runif(2, -8, 8)); gr <- seq(ab[1], ab[2], length.out = 400)
  wg <- n_tr * plogis(gr) * (1 - plogis(gr))
  if (min(wg) < min(wg[1], wg[400]) - 1e-12) bad <- bad + 1
}
cat("intervals where the min is NOT at an endpoint:", bad, "of 2000\n")

## ---- the certified set, and eta ranges over its projection -----------------

negll <- function(b) { s <- 0
  for (j in 1:J) { e <- eta_j(j, b[idx(j)])
    s <- s - sum(yl[[j]] * e - n_tr * log1p(exp(e))) }; s }
f    <- function(b) negll(b) + 0.5 * drop(crossprod(b - mu_b, Lam_b %*% (b - mu_b)))
gr_f <- function(b) { g <- drop(Lam_b %*% (b - mu_b))
  for (j in 1:J) { e <- eta_j(j, b[idx(j)])
    g[idx(j)] <- g[idx(j)] - drop(t(Zl[[j]]) %*% (yl[[j]] - n_tr * plogis(e))) }; g }
He_f <- function(b) { H <- Lam_b
  for (j in 1:J) H[idx(j), idx(j)] <- H[idx(j), idx(j)] + G_j_mat(j, b[idx(j)]); H }
newton <- function(b, tgt = gr_f) { for (i in 1:300) {
  s <- solve(He_f(b), tgt(b)); b <- b - s; if (max(abs(s)) < 1e-11) break }; b }

b_dag <- newton(mu_b); f_dag <- f(b_dag)
Xi <- function(b) f(b) - f_dag
r_lev <- 6

supp <- function(cvec) {                 # max c'beta over {Xi <= r_lev}
  gt <- function(t) { b <- newton(b_dag, function(x) gr_f(x) - t * cvec)
                      list(b = b, xi = Xi(b)) }
  hi <- 1; while (gt(hi)$xi < r_lev && hi < 1e8) hi <- hi * 2
  t <- uniroot(function(t) gt(t)$xi - r_lev, c(1e-10, hi), tol = 1e-11)$root
  drop(crossprod(cvec, gt(t)$b))
}
emb <- function(j, u) { v <- numeric(n_dim); v[idx(j)] <- u; v }
wof <- function(e) n_tr * plogis(e) * (1 - plogis(e))

## ---- assemble the two floors -----------------------------------------------

fl <- lapply(1:J, function(j) {
  Z <- Zl[[j]]
  hi <-  sapply(1:nobs, function(i) supp(emb(j,  Z[i, ])))
  lo <- -sapply(1:nobs, function(i) supp(emb(j, -Z[i, ])))
  wb <- pmin(wof(lo), wof(hi))            # fact (ii): endpoints suffice
  list(eta_lo = lo, eta_hi = hi, wbar = wb,
       G1 = min(wb) * crossprod(Z), G2 = t(Z) %*% (wb * Z))
})

## points inside B_r: walk out along a random ray to a random fraction of t_r
samp_Br <- function(m) t(sapply(1:m, function(k) {
  u <- rnorm(n_dim); u <- u / sqrt(sum(u^2))
  hi <- 1; while (Xi(b_dag + hi * u) < r_lev) hi <- hi * 2
  tr <- uniroot(function(t) Xi(b_dag + t * u) - r_lev, c(0, hi), tol = 1e-10)$root
  b_dag + runif(1)^(1 / n_dim) * tr * u
}))
S <- samp_Br(600)

cat("\n=== Validity and tightness of the matrix floors (r =", r_lev, ") ===\n")
om <- function(G) max(eigen(solve(Pb + G, Pb), only.values = TRUE)$values)
res <- do.call(rbind, lapply(1:J, function(j) {
  v1 <- v2 <- 0
  for (k in 1:nrow(S)) {
    Gk <- G_j_mat(j, S[k, idx(j)])
    if (min(eigen(Gk - fl[[j]]$G1, TRUE, only.values = TRUE)$values) < -1e-8) v1 <- v1 + 1
    if (min(eigen(Gk - fl[[j]]$G2, TRUE, only.values = TRUE)$values) < -1e-8) v2 <- v2 + 1
  }
  e2 <- eigen(fl[[j]]$G2, symmetric = TRUE, only.values = TRUE)$values
  ez <- eigen(crossprod(Zl[[j]]), symmetric = TRUE, only.values = TRUE)$values
  data.frame(j = j, viol_G1 = v1, viol_G2 = v2,
             wbar_min = min(fl[[j]]$wbar), wbar_max = max(fl[[j]]$wbar),
             G2_lmin = e2[2], G2_lmax = e2[1],
             cond_G2 = e2[1] / e2[2], cond_ZtZ = ez[1] / ez[2],
             om_G1 = om(fl[[j]]$G1), om_G2 = om(fl[[j]]$G2),
             G2_dom_G1 = min(eigen(fl[[j]]$G2 - fl[[j]]$G1, TRUE,
                                   only.values = TRUE)$values) >= -1e-10)
}))
print(res, row.names = FALSE, digits = 3)

cat("\nPSD violations over", nrow(S), "points in B_r (must be 0):",
    sum(res$viol_G1), "and", sum(res$viol_G2), "\n")
cat("Gamma^2 >= Gamma^1 for every group:", all(res$G2_dom_G1), "\n")
cat(sprintf("omega_max:  scalar floor %.4f   per-observation floor %.4f\n",
            max(res$om_G1), max(res$om_G2)))
cat("\nWithin-group weight spread wbar_max/wbar_min:",
    paste(format(res$wbar_max / res$wbar_min, digits = 3), collapse = " "), "\n")
cat("That spread is exactly what Gamma^1 throws away by using one worst weight.\n")
cat("\nThe floor can be MORE anisotropic than the design: cond(Gamma^2) vs",
    "\ncond(Z'Z) above, since observations differ in how far their own eta can",
    "\nrun inside B_r.  omega_max is set by the worst DIRECTION in the worst",
    "\ngroup, not by a worst group alone.\n")

## ---- Which of these problems are actually convex? --------------------------
##
## Three different optimisations appear, and they are NOT all convex:
##
##  (A) support function:  min  -c'beta   s.t.  Xi(beta) <= r
##      linear objective, convex constraint  =>  a genuine convex program.
##      Slater holds (Xi(b_dag) = 0 < r), so KKT is necessary AND sufficient:
##      grad Xi(beta) = c/t, Xi = r.  That is exactly what supp() solves, so its
##      answer is the GLOBAL max.  Checked against ray-shooting below.
##
##  (B) the weight minimum:  min_{beta in B} w(eta(beta))
##      w is log-concave, so this is CONCAVE-type minimisation -- not convex,
##      and in general hard.  We never solve it as posed.  Because w depends on
##      beta only through the scalar eta, and the image of a convex set under a
##      linear map is an interval, it collapses to min{w(eta^-), w(eta^+)} with
##      the endpoints from (A).  The reduction is what makes it easy, not
##      convexity of the objective.
##
##  (C) the directional bound:  min_{beta in B} sum_i w_i(beta) (z_i'v)^2
##      a SUM of log-concave functions: neither convex nor quasi-concave.
##      Genuinely global.  This is why option (3) is not offered.

cat("\n=== Convexity audit ===\n")

## (A1) is Xi convex?  Hessian must dominate Lam_b > 0 everywhere.
lmin_Lam <- min(eigen(Lam_b, symmetric = TRUE, only.values = TRUE)$values)
lmin_H <- min(sapply(1:400, function(k) {
  b <- b_dag + 3 * rnorm(n_dim)
  min(eigen(He_f(b), symmetric = TRUE, only.values = TRUE)$values)
}))
cat(sprintf("lambda_min(Lam_b) = %.4f;  min over 400 random beta of",
            lmin_Lam),
    sprintf("lambda_min(Hess Xi) = %.4f\n", lmin_H))
cat("  Xi strictly convex:", lmin_H >= lmin_Lam - 1e-8, "\n")

## (A2) is the support-function answer the global max?  Shoot rays and compare.
gl <- sapply(1:J, function(j) {
  cv  <- emb(j, Zl[[j]][1, ])
  lag <- supp(cv)
  best <- max(sapply(1:4000, function(k) {
    u <- rnorm(n_dim); u <- u / sqrt(sum(u^2))
    hi <- 1; while (Xi(b_dag + hi * u) < r_lev) hi <- hi * 2
    tr <- uniroot(function(t) Xi(b_dag + t * u) - r_lev, c(0, hi), tol = 1e-9)$root
    drop(crossprod(cv, b_dag + tr * u))
  }))
  c(lagrangian = lag, ray_best = best)
})
cat("\nsupport function: Lagrangian/KKT vs best of 4000 boundary rays\n")
print(round(t(gl), 4))
cat("  KKT value >= every sampled boundary point:",
    all(gl["lagrangian", ] >= gl["ray_best", ] - 1e-6), "\n")

## (B) the weight objective is NOT convex -- confirm it fails the chord test
bad_cvx <- 0
for (k in 1:2000) {
  a <- runif(2, -6, 6); lam <- runif(1)
  m <- lam * a[1] + (1 - lam) * a[2]
  if (wof(m) > lam * wof(a[1]) + (1 - lam) * wof(a[2]) + 1e-12) bad_cvx <- bad_cvx + 1
}
cat(sprintf("\nw(eta) violates convexity on %d of 2000 random chords",
            bad_cvx), "(expected: it is log-CONCAVE)\n")

## (C) is a sum of two logit weights quasi-concave?  A counterexample suffices.
## Two humps at eta = +-2.5 with a dip between them: on the interval joining the
## humps the minimum is INTERIOR, so the endpoint rule of (B) fails for sums.
s2 <- function(e) wof(e - 2.5) + wof(e + 2.5)
gr <- seq(-2.5, 2.5, length.out = 2001); v2 <- s2(gr)
cat(sprintf("\nsum of two shifted weights on [-2.5, 2.5]:"),
    sprintf("min %.4f at eta = %.2f;  endpoints %.4f, %.4f\n",
            min(v2), gr[which.min(v2)], v2[1], v2[2001]))
cat("  minimum is interior:", !(which.min(v2) %in% c(1, 2001)),
    "-> quasi-concavity FAILS for sums, so (C) cannot be solved by\n",
    "  scanning the boundary and is a genuine global optimisation.\n")

## (D) does convexity of -loglik survive other links?
cat("\nCurvature of the per-observation negative log-likelihood in eta:\n")
ee <- seq(-6, 6, length.out = 4001); h <- ee[2] - ee[1]
fam <- list(
  logit   = function(e) log1p(exp(e)) - 0.3 * e,
  probit  = function(e) -(0.3 * pnorm(e, log.p = TRUE) +
                          0.7 * pnorm(-e, log.p = TRUE)),
  cloglog = function(e) -(0.3 * log(-expm1(-exp(e))) - 0.7 * exp(e)),
  pois_log = function(e) exp(e) - 0.3 * e,
  gamma_log = function(e) e + 1.4 * exp(-e))
for (nm in names(fam)) {
  d2 <- diff(fam[[nm]](ee), differences = 2) / h^2
  cat(sprintf("  %-10s min second derivative = %+.4e   convex: %s\n",
              nm, min(d2), min(d2) >= -1e-6))
}

## ---- The feasible set for beta_j depends on the other groups ---------------
##
## B_j := proj_j(B_r) is NOT a group-local object.  Xi couples the groups through
## Lam_b, which is not block diagonal (integrating gamma out couples them), so
##
##   B_j = { beta_j : Xi_j^prof(beta_j) <= r },
##   Xi_j^prof(beta_j) := min_{beta_{-j}} Xi(beta_j, beta_{-j} )   [partial min]
##
## convex, since partial minimisation of a jointly convex function is convex.
## Maximising z'beta_j over B_j is therefore a constrained optimisation over the
## FULL vector, with beta_{-j} free to move so as to pay as little Xi-budget as
## possible.  supp() above already does this: it solves grad Xi(beta) = t c in
## R^n, so the other groups accommodate optimally.
##
## The tempting shortcut -- freeze beta_{-j} at the mode and use the CONDITIONAL
## deficiency Xi(., b_dag_{-j}) -- gives a SMALLER set, since profile <= cond,
## hence an over-optimistic floor.  Same unsafe direction as simulation.

cat("\n=== Profile vs conditional: what the other groups contribute ===\n")

cond_supp <- function(j, u) {          # max u'beta_j with beta_{-j} frozen at mode
  gj <- function(bj) { b <- b_dag; b[idx(j)] <- bj
    drop(gr_f(b)[idx(j)]) }
  Hj <- function(bj) { b <- b_dag; b[idx(j)] <- bj
    He_f(b)[idx(j), idx(j)] }
  xj <- function(bj) { b <- b_dag; b[idx(j)] <- bj; Xi(b) }
  nwt <- function(bj, t) { for (k in 1:300) {
    s <- solve(Hj(bj), gj(bj) - t * u); bj <- bj - s
    if (max(abs(s)) < 1e-11) break }; bj }
  hi <- 1; while (xj(nwt(b_dag[idx(j)], hi)) < r_lev && hi < 1e8) hi <- hi * 2
  t <- uniroot(function(t) xj(nwt(b_dag[idx(j)], t)) - r_lev,
               c(1e-10, hi), tol = 1e-11)$root
  drop(crossprod(u, nwt(b_dag[idx(j)], t)))
}

Lbi <- solve(Lam_b)
prof <- do.call(rbind, lapply(1:J, function(j) {
  Z <- Zl[[j]]
  ## eta half-extent about the mode, profile (true) vs conditional (frozen)
  ph <- sapply(1:nobs, function(i) supp(emb(j, Z[i, ]))) -
        drop(Z %*% b_dag[idx(j)])
  ch <- sapply(1:nobs, function(i) cond_supp(j, Z[i, ])) -
        drop(Z %*% b_dag[idx(j)])
  ## ellipsoid factors: marginal block of the inverse vs inverse of the block
  mrg <- mean(sapply(1:nobs, function(i)
    sqrt(drop(crossprod(Z[i, ], Lbi[idx(j), idx(j)] %*% Z[i, ])))))
  cnd <- mean(sapply(1:nobs, function(i)
    sqrt(drop(crossprod(Z[i, ], solve(Lam_b[idx(j), idx(j)], Z[i, ]))))))
  wb_c <- pmin(wof(drop(Z %*% b_dag[idx(j)]) - ch),
               wof(drop(Z %*% b_dag[idx(j)]) + ch))
  data.frame(j = j, eta_prof = mean(ph), eta_cond = mean(ch),
             ratio = mean(ph / ch),
             ell_marg = mrg, ell_cond = cnd, ell_ratio = mrg / cnd,
             om_prof = om(fl[[j]]$G2),
             om_cond_UNSAFE = om(t(Z) %*% (wb_c * Z)))
}))
print(prof, row.names = FALSE, digits = 4)
cat("\nProfile extents exceed conditional ones by a factor of",
    format(mean(prof$ratio), digits = 3),
    "on average,\nso freezing the other groups would understate the excursion",
    "and overstate\nthe floor: omega_max", format(max(prof$om_cond_UNSAFE), digits = 4),
    "(unsafe) vs", format(max(prof$om_prof), digits = 4), "(correct).\n")
cat("\nThe ellipsoid relaxation (F) already gets this right: it uses the",
    "\nMARGINAL block (Lam_b^{-1})_jj, not ((Lam_b)_jj)^{-1}.  Ratio of the two:",
    format(mean(prof$ell_ratio), digits = 4), "\n")

## ---- Can the minima be found by simulating from B_r instead? ----------------
##
## Two effects pull opposite ways:
##   (a) ASSEMBLY SLACK (safe).  Gamma^2 uses each observation's own minimum,
##       and those are not simultaneously attainable, so Gamma^2 sits strictly
##       below the tightest valid floor.
##   (b) SIMULATION ERROR (unsafe).  min over draws >= true min, so a simulated
##       weight floor is too LARGE, making Gamma too large -- i.e. not a lower
##       bound at all.
## The question is whether (a) absorbs (b).  Tested adversarially: the exact
## extremal points from the support function are the places where a simulated
## floor should fail, so those are the test points.

cat("\n=== Simulated minima vs exact support function ===\n")

extremal <- function(cvec) {              # the argmax point, not just the value
  gt <- function(t) { b <- newton(b_dag, function(x) gr_f(x) - t * cvec)
                      list(b = b, xi = Xi(b)) }
  hi <- 1; while (gt(hi)$xi < r_lev && hi < 1e8) hi <- hi * 2
  t <- uniroot(function(t) gt(t)$xi - r_lev, c(1e-10, hi), tol = 1e-11)$root
  gt(t)$b
}
adv <- do.call(rbind, lapply(1:J, function(j)
  do.call(rbind, lapply(1:nobs, function(i)
    rbind(extremal(emb(j,  Zl[[j]][i, ])), extremal(emb(j, -Zl[[j]][i, ])))))))

sim_report <- function(M) {
  D <- samp_Br(M)
  out <- do.call(rbind, lapply(1:J, function(j) {
    Z <- Zl[[j]]; E <- D[, idx(j), drop = FALSE] %*% t(Z)
    hi_s <- apply(E, 2, max); lo_s <- apply(E, 2, min)
    wb_s <- pmin(wof(lo_s), wof(hi_s))
    G2s  <- t(Z) %*% (wb_s * Z)
    rng_true <- fl[[j]]$eta_hi - fl[[j]]$eta_lo
    v <- sum(sapply(1:nrow(adv), function(k)
      min(eigen(G_j_mat(j, adv[k, idx(j)]) - G2s, TRUE,
                only.values = TRUE)$values) < -1e-8))
    data.frame(j = j, range_recovered = mean((hi_s - lo_s) / rng_true),
               wbar_ratio = mean(wb_s / fl[[j]]$wbar),
               viol_adv = v, om_sim = om(G2s), om_exact = om(fl[[j]]$G2))
  }))
  cat(sprintf("\nM = %d draws\n", M)); print(out, row.names = FALSE, digits = 3)
  cat(sprintf("  total PSD violations at extremal points: %d of %d;",
              sum(out$viol_adv), J * nrow(adv)),
      sprintf(" omega_max sim %.4f vs exact %.4f\n",
              max(out$om_sim), max(out$om_exact)))
  invisible(out)
}
for (M in c(200L, 2000L, 20000L)) sim_report(M)

## ---- how the shortfall scales with dimension -------------------------------
##
## Cheap, model-free version of the same effect: uniform draws in the unit ball
## of R^n, recovering the maximum of a FIXED linear functional (true max = 1).
## This is the geometry the simulation route is up against, and it is the reason
## the shortfall above cannot be fixed by more draws.

cat("\n=== Why more draws will not fix it: max of a fixed direction ===\n")
rec <- do.call(rbind, lapply(c(2, 5, 10, 20, 50, 100), function(n) {
  row <- sapply(c(1e3, 1e5), function(M) {
    M <- as.integer(M)
    ## uniform in the n-ball: first coordinate has density prop to (1-x^2)^((n-1)/2)
    x <- rnorm(M); r <- runif(M)^(1 / n)
    g <- matrix(rnorm(M * n), M, n)
    max(r * g[, 1] / sqrt(rowSums(g^2)))
  })
  data.frame(n = n, M_1e3 = row[1], M_1e5 = row[2])
}))
print(rec, row.names = FALSE, digits = 3)
cat("\nFraction of the true extent (= 1) recovered.  The cost of an extra digit",
    "\ngrows like M^(2/(n+1)), so in the dimensions that matter (n = J p_re)",
    "\nsimulation under-explores exactly the directions the floor depends on.\n")
