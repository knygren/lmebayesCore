## =============================================================================
## Exact certified curvature floor for a logit GLMM:  J = 10 groups,
## p_re = 3 random-effect dimensions per group  (n = J * p_re = 30),
## in the WEAK-LIKELIHOOD / HEAVY-SHRINKAGE regime.
##
## Everything below is a Newton solve or a closed form.  There is no simulation:
## the minima are attained exactly, not estimated.
##
##   wbar_ij = min { w(z_ij' beta_j) : beta in Btilde_r },  Btilde_r = {Xi <= r}
##
## is computed by the two-step reduction
##
##   min_{beta in B} w(c'beta) = min_{eta in [eta-,eta+]} w(eta)
##                             = min{ w(eta-), w(eta+) }                    (EP)
##
## with [eta-, eta+] = c' Btilde_r the image of the JOINT set (an interval,
## since Btilde_r is convex compact).  Endpoints come from the KKT system
##
##   grad_j Xi = t z_ij ,      grad_k Xi = 0  (k != j)
##
## so the other groups sit at their CONDITIONAL mode given beta_j -- the profile
## minimiser.  Weak data is exactly the regime where that matters most.
##
## KEY IDENTITY exploited throughout.  In the two-block Gibbs sampler,
##   E[beta_j | gamma] = (G_j + P_b)^{-1} (G_j betahat_j + P_b gamma),
## so (G_j + P_b)^{-1} P_b is BOTH the shrinkage weight on the prior AND the
## coefficient governing how the beta-block follows the gamma-block.  Hence
##   omega_j := lambda_max( (Gamma_j + P_b)^{-1} P_b )
## is simultaneously the certified contraction factor and a certified UPPER
## bound on shrinkage.  "Lots of shrinkage" and "omega near 1" are the same
## statement, so the regime the user asked for is intrinsically the hard one.
## =============================================================================

## ---- fixed across every case: design, truth, within-group prior -------------
set.seed(2026)
J <- 10L; p <- 3L; nobs <- 10L
n_dim <- J * p
idx <- function(j) ((j - 1L) * p + 1L):(j * p)

Zl <- lapply(1:J, function(j)
  cbind(1, matrix(rnorm(nobs * (p - 1L)), nobs, p - 1L)))

Rb     <- matrix(c(1, .30, -.20, .30, 1, .25, -.20, .25, 1), p, p)
sb     <- c(.75, .55, .45)
Sig_b  <- diag(sb) %*% Rb %*% diag(sb)
Lam_b1 <- solve(Sig_b)                       # P_b: within-group prior precision
mu0    <- c(-0.4, 0.3, -0.2)
gam_tr <- mu0 + 0.25 * rnorm(p)
b_tr   <- as.vector(vapply(1:J, function(j)
  gam_tr + drop(t(chol(Sig_b)) %*% rnorm(p)), numeric(p)))

## stable pieces
l1pe <- function(e) ifelse(e > 30, e, log1p(exp(e)))
## w(eta) = n p(1-p) = n/(4 cosh^2(eta/2)); the log form avoids cosh overflow,
## which matters here because weak data pushes |eta| into the tens.
wof_n <- function(e, n) n * exp(-abs(e) - 2 * l1pe(-abs(e)))

ePb <- eigen(Lam_b1, symmetric = TRUE)
Pmh <- ePb$vectors %*% diag(1 / sqrt(ePb$values)) %*% t(ePb$vectors)
## omega = lambda_max((G+P_b)^{-1}P_b); the matrix is NOT symmetric, so use
## G = A M A with A = P_b^{1/2}, giving omega = 1/(1 + lambda_min(M)).
omega_of <- function(Gm) 1 / (1 + min(eigen(Pmh %*% Gm %*% Pmh,
                                            symmetric = TRUE)$values))

ray_bound <- function(r) pgamma(r, n_dim, lower.tail = FALSE) /
                         pgamma(r, n_dim, lower.tail = TRUE)
r_for <- function(eps) uniroot(function(r) ray_bound(r) - eps,
                               c(n_dim, 60 * n_dim), tol = 1e-9)$root

## ---- one case = (trials per observation, population prior precision) --------

make_case <- function(n_tr, lam_g_s, eps = 0.01, seed = 7L) {
  set.seed(seed)
  yl <- lapply(1:J, function(j)
    rbinom(nobs, n_tr, plogis(drop(Zl[[j]] %*% b_tr[idx(j)]))))
  Lam_g <- diag(rep(lam_g_s, p))
  Mcpl  <- Lam_b1 %*% solve(Lam_g + J * Lam_b1) %*% Lam_b1
  Lam_b <- kronecker(diag(J), Lam_b1) - kronecker(matrix(1, J, J), Mcpl)
  mu_b  <- rep(mu0, J)

  eta_j   <- function(j, bj) drop(Zl[[j]] %*% bj)
  negll_j <- function(j, bj) { e <- eta_j(j, bj)
    -sum(yl[[j]] * e - n_tr * l1pe(e)) }
  quadp <- function(b) 0.5 * drop(crossprod(b - mu_b, Lam_b %*% (b - mu_b)))
  f     <- function(b) sum(vapply(1:J, function(j) negll_j(j, b[idx(j)]), 0)) +
                       quadp(b)
  gr_f <- function(b) { g <- drop(Lam_b %*% (b - mu_b))
    for (j in 1:J) { e <- eta_j(j, b[idx(j)])
      g[idx(j)] <- g[idx(j)] - drop(crossprod(Zl[[j]], yl[[j]] - n_tr * plogis(e))) }
    g }
  Gmat <- function(j, bj) { e <- eta_j(j, bj)
    crossprod(Zl[[j]], (n_tr * plogis(e) * (1 - plogis(e))) * Zl[[j]]) }
  He_f <- function(b) { H <- Lam_b
    for (j in 1:J) H[idx(j), idx(j)] <- H[idx(j), idx(j)] + Gmat(j, b[idx(j)]); H }
  newton <- function(b0, rhs = function(x) 0) { b <- b0
    for (k in 1:300) { s <- solve(He_f(b), gr_f(b) - rhs(b))
      st <- 1; while (st > 1e-10 && max(abs(s * st)) > 5) st <- st / 2
      b <- b - st * s; if (max(abs(s)) < 1e-10) break }
    b }

  b_dag <- newton(mu_b); f_dag <- f(b_dag)
  Xi <- function(b) f(b) - f_dag
  r_lev <- r_for(eps)

  emb  <- function(j, u) { v <- numeric(n_dim); v[idx(j)] <- u; v }
  wof  <- function(e) wof_n(e, n_tr)
  pathb <- function(t, cv, start) newton(start, function(x) t * cv)
  supp <- function(cv, r = r_lev, start = b_dag) {
    hi <- 1; bh <- pathb(hi, cv, start)
    while (Xi(bh) < r && hi < 1e12) { hi <- hi * 4; bh <- pathb(hi, cv, bh) }
    tt <- uniroot(function(t) Xi(pathb(t, cv, start)) - r, c(0, hi), tol = 1e-10)$root
    b <- pathb(tt, cv, start)
    list(val = drop(crossprod(cv, b)), beta = b, t = tt) }
  supp_frozen <- function(j, u, r = r_lev) {
    gj <- function(bj) { b <- b_dag; b[idx(j)] <- bj; gr_f(b)[idx(j)] }
    Hj <- function(bj) { b <- b_dag; b[idx(j)] <- bj; He_f(b)[idx(j), idx(j)] }
    xj <- function(bj) { b <- b_dag; b[idx(j)] <- bj; Xi(b) }
    nw <- function(t) { bj <- b_dag[idx(j)]
      for (k in 1:300) { s <- solve(Hj(bj), gj(bj) - t * u); bj <- bj - s
        if (max(abs(s)) < 1e-10) break }; bj }
    hi <- 1; while (xj(nw(hi)) < r && hi < 1e12) hi <- hi * 4
    tt <- uniroot(function(t) xj(nw(t)) - r, c(0, hi), tol = 1e-10)$root
    drop(crossprod(u, nw(tt))) }

  environment()
}

floor_of <- function(cs, r = cs$r_lev, frozen = FALSE) lapply(1:J, function(j) {
  Z <- Zl[[j]]
  if (frozen) {
    hi <-  vapply(1:nobs, function(i) cs$supp_frozen(j,  Z[i, ], r), 0)
    lo <- -vapply(1:nobs, function(i) cs$supp_frozen(j, -Z[i, ], r), 0)
  } else {
    hi <-  vapply(1:nobs, function(i) cs$supp(cs$emb(j,  Z[i, ]), r)$val, 0)
    lo <- -vapply(1:nobs, function(i) cs$supp(cs$emb(j, -Z[i, ]), r)$val, 0)
  }
  es <- pmax(abs(lo), abs(hi)); wb <- cs$wof(es)
  list(lo = lo, hi = hi, es = es, wbar = wb, Gam = crossprod(Z, wb * Z))
})

## ---- 1. the featured case: weak data, heavy shrinkage -----------------------

n_tr0 <- 3L; lg0 <- 0.05
cs <- make_case(n_tr0, lg0)
cat(sprintf("logit GLMM: J = %d, p_re = %d, n = %d; %d obs/group, %d trials/obs\n",
            J, p, n_dim, nobs, n_tr0))
cat(sprintf("population prior precision Lam_g = %.3f I  (weak)\n", lg0))
cat(sprintf("mode: ||grad|| = %.1e;  lambda_min(Hess) = %.4f\n",
            max(abs(cs$gr_f(cs$b_dag))), min(eigen(cs$He_f(cs$b_dag),
                                                   symmetric = TRUE)$values)))
shr <- vapply(1:J, function(j) omega_of(cs$Gmat(j, cs$b_dag[idx(j)])), 0)
cat(sprintf("ACTUAL shrinkage at the mode, lambda_max((G_j+P_b)^{-1}P_b):\n  %s\n",
            paste(format(shr, digits = 3), collapse = " ")))
cat(sprintf("  mean %.3f, max %.3f -- the posterior mean of beta_j is mostly\n",
            mean(shr), max(shr)),
    " the prior/gamma, which is what 'heavy shrinkage' means.\n")
cat(sprintf("\nmass budget: escape <= 0.01  =>  r = %.3f  (n = %d)\n",
            cs$r_lev, n_dim))

## ---- 2. the path: how the other groups buy excursion ------------------------

j0 <- 1L; i0 <- 1L
c0 <- cs$emb(j0, Zl[[j0]][i0, ]); res0 <- cs$supp(c0)
cat("\n=== The path of beta(t): where the excursion comes from ===\n")
cat(sprintf("group %d, obs %d; eta at the mode = %+.4f\n",
            j0, i0, drop(crossprod(c0, cs$b_dag))))
bp <- cs$b_dag
tabp <- do.call(rbind, lapply(seq(0, res0$t, length.out = 7L), function(t) {
  b <- cs$pathb(t, c0, bp); bp <<- b; d <- b - cs$b_dag
  Lj  <- cs$negll_j(j0, b[idx(j0)]) - cs$negll_j(j0, cs$b_dag[idx(j0)])
  Lmj <- sum(vapply(setdiff(1:J, j0), function(k)
    cs$negll_j(k, b[idx(k)]) - cs$negll_j(k, cs$b_dag[idx(k)]), 0))
  data.frame(t = t, eta = drop(crossprod(c0, b)), Xi = cs$Xi(b),
             d_own = sqrt(sum(d[idx(j0)]^2)), d_others = sqrt(sum(d[-idx(j0)]^2)),
             L_j = Lj, L_others = Lmj, Q = cs$quadp(b) - cs$quadp(cs$b_dag))
}))
print(tabp, row.names = FALSE, digits = 4)
cat(sprintf("identity  max |Xi - (L_j + L_others + Q)| = %.1e\n",
            max(abs(tabp$Xi - (tabp$L_j + tabp$L_others + tabp$Q)))))
cat(sprintf("the others travel %.2f x as far as group %d itself\n",
            tabp$d_others[nrow(tabp)] / tabp$d_own[nrow(tabp)], j0))

d <- res0$beta - cs$b_dag; D <- matrix(d, p, J)
Dm <- D[, -j0, drop = FALSE]; com <- rowMeans(Dm)
cat(sprintf("\nshare of the others' displacement in the COMMON direction: %.3f\n",
            sum(com^2) * (J - 1) / sum(Dm^2)))
cat(sprintf("cos angle(own, others' mean) = %+.4f\n",
            sum(D[, j0] * com) / sqrt(sum(D[, j0]^2) * sum(com^2))))
cat("With weak data the whole cloud drifts together: a common shift is charged\n",
    "only by Lam_g, which is nearly nothing, so the certified set is huge.\n")

## ---- 3. free vs frozen, now that the coupling is not negligible -------------

cat("\n=== Letting the other groups move vs freezing them at the mode ===\n")
cmp <- do.call(rbind, lapply(1:J, function(j) {
  u <- Zl[[j]][1, ]; cv <- cs$emb(j, u); e0 <- drop(crossprod(cv, cs$b_dag))
  hp <- cs$supp(cv)$val - e0; hc <- cs$supp_frozen(j, u) - e0
  data.frame(group = j, eta_mode = e0, half_free = hp, half_frozen = hc,
             gain = hp / hc, wbar_free = cs$wof(abs(e0) + hp),
             wbar_frozen = cs$wof(abs(e0) + hc))
}))
print(cmp, row.names = FALSE, digits = 4)
cat(sprintf("\nfree/frozen extent: mean %.3f, max %.3f\n",
            mean(cmp$gain), max(cmp$gain)))
cat(sprintf("freezing inflates the weight floor by mean %.3f, max %.3f\n",
            mean(cmp$wbar_frozen / cmp$wbar_free),
            max(cmp$wbar_frozen / cmp$wbar_free)))

## How big is the coupling correction?  H_{j,-j} = -(Mcpl repeated), and
## Mcpl = Lam_b (Lam_g + J Lam_b)^{-1} Lam_b ~ Lam_b / J for weak Lam_g, so the
## Schur correction is O(Lam_b / J): DILUTED BY THE NUMBER OF GROUPS.  With
## J = 10 profile and conditional curvature nearly coincide however weak the
## data are.  The tempting formula S = (Lam_g^{-1} + Sig_b)^{-1} is the H =
## Lam_beta limit; it needs the likelihood to be negligible against Lam_beta,
## which is NOT the case here, and it is reported below only to show the gap.
H <- cs$He_f(cs$b_dag); o <- idx(j0)
S <- H[o, o] - H[o, -o] %*% solve(H[-o, -o], H[-o, o])
cat(sprintf("\nprofile (Schur) vs conditional curvature at the mode:\n"))
cat(sprintf("  lambda_min: %.4f vs %.4f   ratio %.3f\n",
            min(eigen(S, symmetric = TRUE)$values),
            min(eigen(H[o, o], symmetric = TRUE)$values),
            min(eigen(H[o, o], symmetric = TRUE)$values) /
            min(eigen(S, symmetric = TRUE)$values)))
cat(sprintf("  the H = Lam_beta limit (Lam_g^{-1}+Sig_b)^{-1} would give %.4f;\n",
            min(eigen(solve(solve(diag(rep(lg0, p))) + Sig_b),
                      symmetric = TRUE)$values)),
    " it does not apply, because the likelihood is not negligible here.\n")
cat(sprintf("  ||Mcpl|| / ||Lam_b1|| = %.4f  (~ 1/J = %.3f)\n",
            norm(Lam_b1 %*% solve(diag(rep(lg0, p)) + J * Lam_b1) %*% Lam_b1, "F") /
            norm(Lam_b1, "F"), 1 / J))

## ---- 4. the floor, and whether it certifies ---------------------------------

cat("\n=== The floor at the featured (weak) case ===\n")
fl <- floor_of(cs)
adv <- do.call(rbind, lapply(1:J, function(j) t(vapply(1:nobs, function(i)
  cs$supp(cs$emb(j, Zl[[j]][i, ]))$beta, numeric(n_dim)))))
viol <- sum(vapply(seq_len(nrow(adv)), function(k) sum(vapply(1:J, function(j)
  as.integer(min(eigen(cs$Gmat(j, adv[k, idx(j)]) - fl[[j]]$Gam,
                       symmetric = TRUE)$values) < -1e-8), 0L)), 0L))
sm <- do.call(rbind, lapply(1:J, function(j)
  data.frame(group = j, eta_star_max = max(fl[[j]]$es),
             wbar_min = min(fl[[j]]$wbar),
             shrink_at_mode = shr[j], omega_certified = omega_of(fl[[j]]$Gam))))
print(sm, row.names = FALSE, digits = 4)
cat(sprintf("\nPSD violations at %d extremal points: %d\n", nrow(adv), viol))
cat(sprintf("omega_max certified = %.6f   (actual shrinkage at the mode: %.4f)\n",
            max(sm$omega_certified), max(shr)))
cat("omega_certified >= shrink_at_mode always, since Gamma_j <= G_j(bdag_j).\n")

## ---- 5. sweep 1: weakening the likelihood -----------------------------------

cat("\n=== Sweep 1: what weakening the data does (Lam_g fixed at 0.05) ===\n")
sw1 <- do.call(rbind, lapply(c(25L, 10L, 3L, 1L), function(nt) {
  ci <- make_case(nt, lg0); fi <- floor_of(ci)
  sh <- vapply(1:J, function(j) omega_of(ci$Gmat(j, ci$b_dag[idx(j)])), 0)
  data.frame(trials = nt,
             lmin_Hess = min(eigen(ci$He_f(ci$b_dag), symmetric = TRUE)$values),
             shrink_max = max(sh),
             eta_star_max = max(vapply(fi, function(z) max(z$es), 0)),
             wbar_min = min(vapply(fi, function(z) min(z$wbar), 0)),
             omega_max = max(vapply(fi, function(z) omega_of(z$Gam), 0)))
}))
print(sw1, row.names = FALSE, digits = 4)
cat("\nShrinkage and omega_max move together, as they must: they are the same\n",
    "matrix. Weak data buys heavy shrinkage AND a vacuous certificate at once.\n")

## ---- 6. sweep 2: what actually rescues it -----------------------------------

cat("\n=== Sweep 2: tightening the POPULATION prior (the hypothesis to kill) ===\n")
sw2 <- do.call(rbind, lapply(c(0.05, 2, 50), function(lg) {
  ci <- make_case(n_tr0, lg); fi <- floor_of(ci)
  sh <- vapply(1:J, function(j) omega_of(ci$Gmat(j, ci$b_dag[idx(j)])), 0)
  data.frame(lam_g = lg,
             marg_curv = min(eigen(solve(solve(diag(rep(lg, p))) + Sig_b),
                                   symmetric = TRUE)$values),
             shrink_max = max(sh),
             eta_star_max = max(vapply(fi, function(z) max(z$es), 0)),
             omega_max = max(vapply(fi, function(z) omega_of(z$Gam), 0)))
}))
print(sw2, row.names = FALSE, digits = 4)
cat("\nA 1000-fold increase in Lam_g moves eta_star by ~4% and leaves omega_max\n",
    "at 1.  The population prior is NOT the binding constraint: its influence\n",
    "reaches beta_j only through Mcpl ~ Lam_b / J, which J = 10 has diluted.\n")

## ---- 7. what IS binding: the level r, and hence the dimension ---------------
## Decouple r from the mass budget and sweep it directly.

cat("\n=== Sweep 3: the level r itself, at 3 trials/obs ===\n")
sw3 <- do.call(rbind, lapply(c(3, 6, 12, 24, cs$r_lev), function(rr) {
  fi <- floor_of(cs, rr)
  data.frame(r = rr, escape_bound = ray_bound(rr),
             eta_star_max = max(vapply(fi, function(z) max(z$es), 0)),
             wbar_min = min(vapply(fi, function(z) min(z$wbar), 0)),
             omega_max = max(vapply(fi, function(z) omega_of(z$Gam), 0)))
}))
print(sw3, row.names = FALSE, digits = 4)

## The mass bound fixes what r must be, and that grows with the dimension.
rq <- vapply(c(6L, 12L, 20L, 30L, 60L), function(nn)
  uniroot(function(r) pgamma(r, nn, lower.tail = FALSE) /
                      pgamma(r, nn, lower.tail = TRUE) - 0.01,
          c(nn, 80 * nn), tol = 1e-9)$root, 0)
print(data.frame(n = c(6L, 12L, 20L, 30L, 60L), r_needed_at_1pct = rq),
      row.names = FALSE, digits = 5)
cat("\nomega_max is useful only for r well below what the mass bound demands at\n",
    "n = 30.  The two requirements pull in opposite directions and the gap is\n",
    "set by the DIMENSION, not by the prior on gamma and not by shrinkage per\n",
    "se: r_needed grows roughly linearly in n = J p_re, while the excursion\n",
    "eta_star ~ sqrt(2r / curvature) grows like sqrt(r).\n")

## ---- 8. how loose is the mass bound? ----------------------------------------
## Proposition 2 assumes only CONVEXITY of Xi.  Its extremal case is a cone,
## Xi growing LINEARLY along rays: substituting v = Xi gives s ~ v and the polar
## Jacobian s^{n-1} ds ~ v^{n-1} dv, hence Gamma(n, r) / gamma(n, r).
## A near-Gaussian posterior has Xi growing QUADRATICALLY along rays, so
## s ~ sqrt(v), s^{n-1} ds ~ v^{n/2-1} dv, and the reference is Gamma(n/2).
## That is a factor of two in the shape parameter -- pure slack for GLM
## posteriors, whose Xi is quadratic-like near the mode.

cat("\n=== How loose is Proposition 2? ===\n")
ray_gauss <- function(r, nn = n_dim) pgamma(r, nn / 2, lower.tail = FALSE) /
                                     pgamma(r, nn / 2, lower.tail = TRUE)
r_gauss <- function(eps, nn = n_dim)
  uniroot(function(r) ray_gauss(r, nn) - eps, c(1, 80 * nn), tol = 1e-9)$root
loose <- data.frame(n = c(6L, 12L, 20L, 30L, 60L),
                    r_convex_only = rq,
                    r_gaussian_ref = vapply(c(6L, 12L, 20L, 30L, 60L),
                                            function(nn) r_gauss(0.01, nn), 0))
loose$ratio <- loose$r_convex_only / loose$r_gaussian_ref
print(loose, row.names = FALSE, digits = 5)
cat(sprintf("\nAt n = %d, Prop 2 demands r = %.2f; the Gaussian reference needs\n",
            n_dim, cs$r_lev),
    sprintf("  only r = %.2f.  At Prop 2's r the Gaussian escape would be %.2e,\n",
            r_gauss(0.01), ray_gauss(cs$r_lev)),
    "  i.e. Prop 2 is spending ~2x the deficiency budget it needs to.\n")

## Does closing that gap rescue the certificate here?
rg <- r_gauss(0.01)
fg <- floor_of(cs, rg)
cat(sprintf("\nomega_max at the (tighter) Gaussian-reference level r = %.2f: %.6f\n",
            rg, max(vapply(fg, function(z) omega_of(z$Gam), 0))))
cat(sprintf("  vs %.6f at Prop 2's r = %.2f\n",
            max(vapply(floor_of(cs), function(z) omega_of(z$Gam), 0)), cs$r_lev))
cat("\nSo the bound IS loose -- by about a factor of two in r -- but tightening\n",
    "it does not rescue this regime.  eta_star ~ sqrt(r) while w decays like\n",
    "exp(-|eta|), so the floor is exp(-C sqrt(r)); halving r buys a factor\n",
    "exp(C sqrt(r)(1 - 1/sqrt2)), which is far too little when omega is\n",
    "already 1 - 4e-5.  The looseness is real and worth fixing; it is not what\n",
    "is killing the n = 30 weak-data case.\n")
cat("\nFixing it rigorously needs MORE than convexity: the cone is a genuine\n",
    "convex worst case, so Gamma(n) cannot be improved without an extra\n",
    "hypothesis (strong convexity, or quadratic domination along rays).  Xi is\n",
    "strongly convex here -- Hess Xi >= Lam_beta > 0 -- so the hypothesis is\n",
    "available; the ray argument under it is not yet written down.\n")
