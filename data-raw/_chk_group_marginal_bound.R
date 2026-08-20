## =============================================================================
## Can we certify GROUP BY GROUP instead of jointly?
##
## Proposition 2 needs only convexity of Xi and the DIMENSION of the space the
## density lives on.  By Prekopa the marginal
##      pitilde_j(beta_j) = int pitilde(beta | y) d beta_{-j}
## is log-concave on R^p, so Prop 2 applies to it in dimension p = 3, not
## n = J p = 30.  A union bound over the J groups costs only log J.  That is the
## prize.  The bill is that the marginal level set is not computable, so we must
## sandwich it by something that is.
##
## THE SANDWICH.  Write, with Xi^prof_j(beta_j) = min_{beta_{-j}} Xi(beta),
##      pitilde_j(beta_j) / pitilde(beta_dag) = exp(-Xi^prof_j(beta_j)) V(beta_j),
##      V(beta_j) = int exp(-[Xi(beta_j, beta_{-j}) - Xi^prof_j(beta_j)]) dbeta_{-j}.
## The bracket is convex in beta_{-j} with minimum 0 and
##      Lam_beta^{-j,-j}  <=  Hess_{-j} Xi  <=  Lam_beta^{-j,-j} + Gbar_{-j},
## Gbar = the max-weight data precision (n_tr/4 for logit).  Hence V lies between
## two Gaussian volumes and, for ANY two arguments,
##      |log V(a) - log V(b)|  <=  kappa_j := 0.5 logdet(I + Lam^{-1} Gbar)_{-j}.
## Two applications (one for the level, one because the marginal mode differs
## from the joint mode) give the CERTIFIED INCLUSION
##
##      { Xi^marg_j <= r }  subset  { Xi^prof_j <= r + 2 kappa_j }.               (M)
##
## The right side is exactly what supp() already computes.  So the whole scheme
## is a change of LEVEL, from r_joint(n) to r_grp(p) + 2 kappa_j.  This script
## asks whether that trade is favourable.
##
## Note the Gaussian sanity check: if Xi is exactly quadratic then Xi^prof_j and
## Xi^marg_j are the SAME Schur complement and kappa_j = 0.  So kappa_j prices
## non-Gaussianity only -- but the bound above prices it crudely.
## =============================================================================

## ---- setup: replicated verbatim from _ex_logit_floor_J10_p3.R ---------------
set.seed(2026)
J <- 10L; p <- 3L; nobs <- 10L
n_dim <- J * p
idx <- function(j) ((j - 1L) * p + 1L):(j * p)

Zl <- lapply(1:J, function(j)
  cbind(1, matrix(rnorm(nobs * (p - 1L)), nobs, p - 1L)))

Rb     <- matrix(c(1, .30, -.20, .30, 1, .25, -.20, .25, 1), p, p)
sb     <- c(.75, .55, .45)
Sig_b  <- diag(sb) %*% Rb %*% diag(sb)
Lam_b1 <- solve(Sig_b)
mu0    <- c(-0.4, 0.3, -0.2)
gam_tr <- mu0 + 0.25 * rnorm(p)
b_tr   <- as.vector(vapply(1:J, function(j)
  gam_tr + drop(t(chol(Sig_b)) %*% rnorm(p)), numeric(p)))

l1pe  <- function(e) ifelse(e > 30, e, log1p(exp(e)))
wof_n <- function(e, n) n * exp(-abs(e) - 2 * l1pe(-abs(e)))

ePb <- eigen(Lam_b1, symmetric = TRUE)
Pmh <- ePb$vectors %*% diag(1 / sqrt(ePb$values)) %*% t(ePb$vectors)
omega_of <- function(Gm) 1 / (1 + min(eigen(Pmh %*% Gm %*% Pmh,
                                            symmetric = TRUE)$values))

make_case <- function(n_tr, lam_g_s, seed = 7L) {
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

  emb  <- function(j, u) { v <- numeric(n_dim); v[idx(j)] <- u; v }
  wof  <- function(e) wof_n(e, n_tr)
  pathb <- function(t, cv, start) newton(start, function(x) t * cv)
  supp <- function(cv, r, start = b_dag) {
    hi <- 1; bh <- pathb(hi, cv, start)
    while (Xi(bh) < r && hi < 1e12) { hi <- hi * 4; bh <- pathb(hi, cv, bh) }
    tt <- uniroot(function(t) Xi(pathb(t, cv, start)) - r, c(0, hi), tol = 1e-10)$root
    b <- pathb(tt, cv, start)
    list(val = drop(crossprod(cv, b)), beta = b, t = tt) }
  environment()
}

floor_at <- function(cs, r) lapply(1:J, function(j) {
  Z <- Zl[[j]]
  hi <-  vapply(1:nobs, function(i) cs$supp(cs$emb(j,  Z[i, ]), r)$val, 0)
  lo <- -vapply(1:nobs, function(i) cs$supp(cs$emb(j, -Z[i, ]), r)$val, 0)
  es <- pmax(abs(lo), abs(hi)); wb <- cs$wof(es)
  list(es = es, wbar = wb, Gam = crossprod(Z, wb * Z))
})

n_tr0 <- 3L; lg0 <- 0.05
cs <- make_case(n_tr0, lg0)

## ---- 1. the prize: what dimension does Proposition 2 actually see? ----------

eps <- 0.01
r_at <- function(d, e) uniroot(function(r) pgamma(r, d, lower.tail = FALSE) /
                                           pgamma(r, d, lower.tail = TRUE) - e,
                               c(1e-8, 200 * d + 200), tol = 1e-10)$root

r_joint <- r_at(n_dim, eps)            # one set in R^30, budget 0.01
r_grp   <- r_at(p,     eps / J)        # J sets in R^3,  budget 0.01/10 each
r_one   <- r_at(1L,    eps / (J * nobs))  # J*nobs sets in R^1, budget 0.01/100

cat("=== the prize: Prop 2 level by the dimension it is applied in ===\n")
pz <- data.frame(
  scheme  = c("joint, R^n", "per group, R^p", "per functional, R^1"),
  dim     = c(n_dim, p, 1L),
  n_sets  = c(1L, J, J * nobs),
  budget  = c(eps, eps / J, eps / (J * nobs)),
  r       = c(r_joint, r_grp, r_one),
  s       = sqrt(2 * c(r_joint, r_grp, r_one)))
print(pz, row.names = FALSE, digits = 4)
cat(sprintf("\nlevel falls %.1fx (joint -> per group) and %.1fx (joint -> per functional);\n",
            r_joint / r_grp, r_joint / r_one),
    sprintf("in radius units s falls %.2fx and %.2fx.  Since the weight floor decays\n",
            sqrt(r_joint / r_grp), sqrt(r_joint / r_one)),
    " like exp(-s/sqrt(kappa)), that is the entire prize.\n")

## ---- 2. the bill: the crude marginalisation constant kappa_j ----------------

Gbar <- lapply(1:J, function(j) crossprod(Zl[[j]], (n_tr0 / 4) * Zl[[j]]))
kap_crude <- vapply(1:J, function(j) {
  o  <- idx(j)
  L  <- cs$Lam_b[-o, -o]
  Gb <- matrix(0, n_dim - p, n_dim - p)
  ks <- setdiff(1:J, j)
  for (a in seq_along(ks)) {
    rg <- ((a - 1L) * p + 1L):(a * p)
    Gb[rg, rg] <- Gbar[[ks[a]]]
  }
  0.5 * as.numeric(determinant(diag(n_dim - p) + solve(L, Gb),
                               logarithm = TRUE)$modulus)
}, 0)

cat("\n=== the bill: crude kappa_j from the uniform Hessian sandwich ===\n")
cat(sprintf("kappa_j (m = %d dims integrated out): min %.2f, max %.2f\n",
            n_dim - p, min(kap_crude), max(kap_crude)))
cat(sprintf("required level R_j = r_grp + 2 kappa_j: min %.2f, max %.2f\n",
            r_grp + 2 * min(kap_crude), r_grp + 2 * max(kap_crude)))
cat(sprintf("compare r_joint = %.2f  ->  crude scheme is %s\n", r_joint,
            ifelse(r_grp + 2 * max(kap_crude) < r_joint, "BETTER", "WORSE")))

## ---- 3. what kappa_j really is: the Laplace value ---------------------------
## log V(beta_j) = (m/2) log 2pi - 0.5 logdet Hess_{-j}(profile point) + o(1).
## Only the VARIATION matters, and it is a determinant RATIO -- which is exactly
## 0 for a Gaussian.  Evaluate it along the actual extremal path.

cat("\n=== what kappa_j really is: determinant variation along the path ===\n")
logdet_mj <- function(b, j) { o <- idx(j)
  as.numeric(determinant(cs$He_f(b)[-o, -o], logarithm = TRUE)$modulus) }
kl <- do.call(rbind, lapply(1:J, function(j) {
  cv <- cs$emb(j, Zl[[j]][1, ])
  be <- cs$supp(cv, r_joint)$beta
  data.frame(group = j, crude = kap_crude[j],
             laplace = 0.5 * abs(logdet_mj(be, j) - logdet_mj(cs$b_dag, j)))
}))
print(kl, row.names = FALSE, digits = 4)
cat(sprintf("\nmean crude %.2f vs mean Laplace %.3f: the sandwich overcharges %.0fx.\n",
            mean(kl$crude), mean(kl$laplace), mean(kl$crude) / mean(kl$laplace)))
cat(sprintf("with the Laplace value, R_j = %.2f + 2(%.3f) = %.2f  vs  r_joint = %.2f\n",
            r_grp, mean(kl$laplace), r_grp + 2 * mean(kl$laplace), r_joint))

## ---- 4. what it would buy: the floor at each level --------------------------

cat("\n=== the floor at each candidate level ===\n")
lv <- c(joint = r_joint,
        grp_crude = r_grp + 2 * max(kap_crude),
        grp_laplace = r_grp + 2 * max(kl$laplace),
        grp_ideal = r_grp,
        one_ideal = r_one)
fo <- do.call(rbind, lapply(names(lv), function(nm) {
  fi <- floor_at(cs, lv[[nm]])
  data.frame(level = nm, r = lv[[nm]], s = sqrt(2 * lv[[nm]]),
             eta_star_max = max(vapply(fi, function(z) max(z$es), 0)),
             wbar_min = min(vapply(fi, function(z) min(z$wbar), 0)),
             omega_max = max(vapply(fi, function(z) omega_of(z$Gam), 0)),
             one_minus_omega = 1 - max(vapply(fi, function(z) omega_of(z$Gam), 0)))
}))
print(fo, row.names = FALSE, digits = 4)
cat("\n1 - omega is the quantity that must be non-negligible for the certificate\n",
    "to say anything.  Read the last column down the table to see whether the\n",
    "dimension reduction is worth the marginalisation constant.\n")

## ---- 5. how much slack is there in kappa? -----------------------------------
## Break-even is R_j = r_joint, i.e. kappa = (r_joint - r_grp)/2.

kap_be <- (r_joint - r_grp) / 2
cat("\n=== how good does the kappa bound have to be? ===\n")
cat(sprintf("break-even kappa (scheme ties the joint bound): %.2f\n", kap_be))
cat(sprintf("crude sandwich gives                           %.2f  (misses by %.0f%%)\n",
            max(kap_crude), 100 * (max(kap_crude) / kap_be - 1)))
cat(sprintf("the Laplace truth is                           %.2f  (%.0fx headroom)\n",
            max(kl$laplace), kap_be / max(kl$laplace)))
kg <- do.call(rbind, lapply(c(0, 0.5, 1, 2, 4, 8, 16), function(kk) {
  fi <- floor_at(cs, r_grp + 2 * kk)
  data.frame(kappa = kk, R = r_grp + 2 * kk, s = sqrt(2 * (r_grp + 2 * kk)),
             one_minus_omega = 1 - max(vapply(fi, function(z) omega_of(z$Gam), 0)),
             gain_vs_joint = (1 - max(vapply(fi, function(z) omega_of(z$Gam), 0))) /
                             4.189e-05)
}))
print(kg, row.names = FALSE, digits = 4)
cat("\nThe payoff is flat in kappa up to about 2 and only then erodes, so a bound\n",
    "anywhere near the truth captures essentially the whole prize.\n")

## ---- 6. does it make the certificate non-vacuous with stronger data? --------

cat("\n=== the same comparison at 25 trials/obs (informative likelihood) ===\n")
cs2 <- make_case(25L, lg0)
st <- do.call(rbind, lapply(list(c("joint", r_joint), c("per group", r_grp + 2)),
                            function(z) {
  fi <- floor_at(cs2, as.numeric(z[2]))
  data.frame(scheme = z[1], r = as.numeric(z[2]), s = sqrt(2 * as.numeric(z[2])),
             eta_star_max = max(vapply(fi, function(q) max(q$es), 0)),
             wbar_min = min(vapply(fi, function(q) min(q$wbar), 0)),
             omega_max = max(vapply(fi, function(q) omega_of(q$Gam), 0)))
}))
print(st, row.names = FALSE, digits = 4)
cat("\nHere omega_max is what decides whether the sampler is certified to\n",
    "contract; the per-group level is what moves it away from 1.\n")

## ---- 7. conditional vs profile vs marginal ---------------------------------
## Three curvatures for beta_j, in increasing variance:
##   conditional  H_jj              others FROZEN at beta_dag_{-j}
##   profile      S_j = Schur       others RELAX to their conditional mode
##   marginal     Hess Xi^marg      others INTEGRATED OUT
## S_j <= H_jj always (Schur), so freezing UNDERSTATES the spread of beta_j.
## For quadratic Xi, marginal == profile exactly.  Beyond that they differ by
## the curvature of the log-volume term 0.5 logdet H_{-j}(betahat_{-j}(beta_j)).

prof_min <- function(cs_, j, bj, start = cs_$b_dag) {
  o <- idx(j); b <- start; b[o] <- bj
  for (k in 1:200) {
    s <- solve(cs_$He_f(b)[-o, -o], cs_$gr_f(b)[-o])
    st <- 1; while (st > 1e-10 && max(abs(s * st)) > 5) st <- st / 2
    b[-o] <- b[-o] - st * s
    if (max(abs(s)) < 1e-10) break
  }
  b
}
ldm <- function(cs_, b, j) { o <- idx(j)
  as.numeric(determinant(cs_$He_f(b)[-o, -o], logarithm = TRUE)$modulus) }

## Laplace marginal deficiency, anchored at the joint mode
xi_marg <- function(cs_, j, bj) { b <- prof_min(cs_, j, bj)
  cs_$Xi(b) + 0.5 * (ldm(cs_, b, j) - ldm(cs_, cs_$b_dag, j)) }
xi_prof <- function(cs_, j, bj) cs_$Xi(prof_min(cs_, j, bj))

hess_fd <- function(fn, x0, h = 1e-3) {
  q <- length(x0); H <- matrix(0, q, q)
  for (a in 1:q) for (b in a:q) {
    ea <- eb <- numeric(q); ea[a] <- h; eb[b] <- h
    H[a, b] <- H[b, a] <- (fn(x0 + ea + eb) - fn(x0 + ea - eb) -
                           fn(x0 - ea + eb) + fn(x0 - ea - eb)) / (4 * h^2)
  }
  H
}

cmp3 <- function(cs_, lab) {
  H <- cs_$He_f(cs_$b_dag)
  out <- do.call(rbind, lapply(1:J, function(j) {
    o  <- idx(j)
    Hc <- H[o, o]
    Sp <- Hc - H[o, -o] %*% solve(H[-o, -o], H[-o, o])
    Hm <- hess_fd(function(x) xi_marg(cs_, j, x), cs_$b_dag[o])
    Hm <- (Hm + t(Hm)) / 2
    z  <- Zl[[j]][1, ]
    sd_of <- function(M) sqrt(drop(crossprod(z, solve(M, z))))
    data.frame(group = j,
               sd_cond = sd_of(Hc), sd_prof = sd_of(Sp), sd_marg = sd_of(Hm),
               prof_over_cond = sd_of(Sp) / sd_of(Hc),
               marg_over_prof = sd_of(Hm) / sd_of(Sp))
  }))
  cat(sprintf("\n=== sd of eta_1j = z' beta_j under the three curvatures (%s) ===\n", lab))
  print(out, row.names = FALSE, digits = 4)
  cat(sprintf("mean profile/conditional = %.4f   mean marginal/profile = %.4f\n",
              mean(out$prof_over_cond), mean(out$marg_over_prof)))
  invisible(out)
}
c3a <- cmp3(cs,  "weak data, 3 trials/obs")
c3b <- cmp3(cs2, "informative data, 25 trials/obs")

cat("\nThe conditional is the anti-conservative one: freezing the other groups\n",
    "shrinks the set for beta_j, RAISES the weight floor, and would certify a\n",
    "contraction the sampler has not earned.  The profile is what supp() gives\n",
    "and is the correct target; the marginal exceeds it only by the curvature of\n",
    "the log-volume term, which is the same non-Gaussian effect priced by kappa.\n")

## Direct check of the level-set consequence: half-extent of eta under each.
cat("\n=== half-extent of eta_1j at r = 11.23, all three (group 1) ===\n")
j <- 1L; u <- Zl[[j]][1, ]; e0 <- drop(crossprod(u, cs$b_dag[idx(j)]))
ext <- function(fn, r) {
  g <- function(t) fn(cs, j, cs$b_dag[idx(j)] + t * u / sum(u^2)) - r
  uniroot(g, c(1e-6, 200), tol = 1e-8, extendInt = "upX")$root
}
cnd <- function(cs_, j, bj) { b <- cs_$b_dag; b[idx(j)] <- bj; cs_$Xi(b) }
cat(sprintf("along the z_1 ray from the mode (eta at mode = %+.3f):\n", e0))
cat(sprintf("  conditional (frozen) : %.3f\n", ext(cnd, r_grp)))
cat(sprintf("  profile              : %.3f\n", ext(xi_prof, r_grp)))
cat(sprintf("  marginal (Laplace)   : %.3f\n", ext(xi_marg, r_grp)))
cat("Ordering conditional <= profile <= marginal is the safety direction:\n",
    "using the conditional would understate how far eta can travel.\n")

## Why the three are so close here, and when they would not be.  The whole
## coupling sits in the off-diagonal blocks of Lam_beta, which are -Mcpl with
##      Mcpl = Lam_b (Lam_g + J Lam_b)^{-1} Lam_b.
## This is monotone decreasing in Lam_g, so it is MAXIMISED at Lam_g -> 0, where
## Mcpl -> Lam_b / J.  The group coupling can therefore never exceed O(1/J), no
## matter how weak the population prior is: with J = 10 the other nine groups
## barely respond to beta_j, and profile, conditional and marginal must agree to
## O(1/J).  Small J (or few groups carrying the information) is where the user's
## concern bites.
cat("\n=== why the three nearly coincide: coupling is O(1/J) ===\n")
dil <- do.call(rbind, lapply(c(1e-6, 0.05, 2, 50), function(lg) {
  Mc <- Lam_b1 %*% solve(diag(rep(lg, p)) + J * Lam_b1) %*% Lam_b1
  data.frame(lam_g = lg, norm_ratio = norm(Mc, "F") / norm(Lam_b1, "F"))
}))
dil$bound_1_over_J <- 1 / J
print(dil, row.names = FALSE, digits = 4)
cat(sprintf("Mcpl / Lam_b is bounded by 1/J = %.3f and attains it as Lam_g -> 0.\n", 1 / J))

## ---- 8. what does the union bound over J actually cost? --------------------
## The per-group scheme already runs at budget eps/J = 0.01/10.  Worth pricing
## separately, because in the TAIL a Gamma(d) quantile grows like log(1/delta)
## for fixed d, so dividing the budget by J costs only about log J in r -- while
## dividing the DIMENSION by J saved a factor of four.  These are the two levers
## and they are not symmetric.

cat("\n=== price of the union bound: budget vs dimension ===\n")
bd <- do.call(rbind, lapply(list(c(3, eps), c(3, eps / J), c(3, eps / J^2),
                                 c(30, eps), c(30, eps / J)), function(z) {
  d <- z[1]; e <- z[2]; r <- r_at(as.integer(d), e)
  data.frame(dim = d, budget = e, r = r, s = sqrt(2 * r))
}))
bd$vs_d3_eps <- bd$r / bd$r[1]
print(bd, row.names = FALSE, digits = 4)
cat(sprintf("\ndividing the budget by J costs  r: %.2f -> %.2f  (+%.2f, log J = %.2f)\n",
            bd$r[1], bd$r[2], bd$r[2] - bd$r[1], log(J)))
cat(sprintf("dividing it by J^2 costs        r: %.2f -> %.2f  (+%.2f more)\n",
            bd$r[2], bd$r[3], bd$r[3] - bd$r[2]))
cat(sprintf("but dropping dim 30 -> 3 saves  r: %.2f -> %.2f  (-%.2f)\n",
            bd$r[4], bd$r[1], bd$r[4] - bd$r[1]))
cat("So the union bound is cheap (additive, ~log J) and the dimension is dear\n",
    "(multiplicative).  Tightening the budget further is affordable.\n")

cat("\n=== effect on the certificate of tightening the total budget ===\n")
tb <- do.call(rbind, lapply(c(eps, eps / J, eps / J^2), function(e) {
  R <- r_at(p, e / J) + 2      # per-group level with a kappa allowance of 1
  fi <- floor_at(cs, R)
  om <- max(vapply(fi, function(z) omega_of(z$Gam), 0))
  data.frame(total_budget = e, per_group = e / J, R = R, s = sqrt(2 * R),
             eta_star_max = max(vapply(fi, function(z) max(z$es), 0)),
             omega_max = om, one_minus_omega = 1 - om)
}))
print(tb, row.names = FALSE, digits = 4)
cat(sprintf("\ntotal budget 0.01 -> %.4f costs only %.1f%% of (1 - omega):\n",
            eps / J, 100 * (1 - tb$one_minus_omega[2] / tb$one_minus_omega[1])),
    "a 10x safer certificate for a small fraction of the strength.\n")

cat("\nSO: the marginal being WIDER is a COST, not a benefit -- at a fixed level a\n",
    "wider set means a lower weight floor.  The benefit of going marginal is\n",
    "entirely the DIMENSION in Proposition 2 (3 instead of 30, r 11.2 instead of\n",
    "44.2).  The accounting is favourable precisely because the width cost is\n",
    "O(1/J) while the level gain is a factor of four.\n")
