# Concrete demonstration: can a "Zellner-like" vector-pwt prior (per-coefficient
# pwt in Prior_Setup()'s vector-pwt formula, Sigma = V0 * outer(s,s) with
# s_i = sqrt((1-pwt_i)/pwt_i), preserving V0's correlation structure exactly)
# produce a genuine interior UB2 minimum that undercuts the endpoint-only
# estimate -- i.e. is this a real, reachable failure mode of a documented
# glmbayes/glmbayesCore feature, not just an adversarially constructed example?
#
# Context: kappa(K) > 2 only means an interior critical point is *not ruled
# out* by Claim 7 part 1/the corrected Remark 5.5.7 bound -- it does not mean
# one necessarily exists. Whether it does depends on the specific face's
# direction (w_i = squared coordinates of v_j in K's eigenbasis) and Delta.
# This script searches over realistic vector-pwt spreads and realistic grid
# faces (mimicking the actual +/- omega per-dimension tangent points used by
# EnvelopeBuild) to find concrete (K, v_j, Delta) triples where the endpoint
# method actually fails, using the exact reduced-form UB2 formula validated
# against the real C++ RSS/UB2 computation in ub2_root_finding_prototype.R.
#
#   Rscript data-raw/vector_pwt_ub2_example.R

devtools::load_all(".", quiet = TRUE)

g_fun     <- function(t, lambda, w) sapply(t, function(tt) sum(w / (lambda + tt)^2))
ub2_red   <- function(t, lambda, w, Delta) (t / 2) * (g_fun(t, lambda, w) - Delta)
hprime    <- function(t, lambda, w) sapply(t, function(tt) sum(w * (lambda - tt) / (lambda + tt)^3))

find_ub2_min_by_roots <- function(lambda, w, Delta, t_lo, t_hi, grid_mult = 60) {
  lam_max <- max(lambda)
  hi_search <- min(t_hi, lam_max * (1 - 1e-9))
  cands <- c(t_lo, t_hi)
  if (hi_search > t_lo) {
    anchors <- sort(unique(c(t_lo, hi_search, lambda[lambda > t_lo & lambda < hi_search])))
    grid <- unique(sort(unlist(lapply(seq_len(max(length(anchors) - 1, 1)), function(i) {
      if (length(anchors) < 2) return(seq(t_lo, hi_search, length.out = grid_mult))
      lo_i <- anchors[i]; hi_i <- anchors[i + 1]
      if (hi_i <= lo_i) return(numeric(0))
      seq(lo_i, hi_i, length.out = grid_mult)
    }))))
    fvals <- sapply(grid, function(t) hprime(t, lambda, w) - Delta)
    sc <- which(fvals[-1] * fvals[-length(fvals)] < 0)
    roots <- vapply(sc, function(i) uniroot(function(t) hprime(t, lambda, w) - Delta,
                                             lower = grid[i], upper = grid[i + 1], tol = 1e-12)$root,
                     numeric(1))
    cands <- c(cands, roots[roots > t_lo & roots < t_hi])
  }
  vals <- sapply(cands, function(t) ub2_red(t, lambda, w, Delta))
  list(t_min = cands[which.min(vals)], ub2_min = min(vals))
}

## ============================================================================
## Step 1: a real, modestly-correlated 2-predictor design (n=25) and its OLS
## quantities, exactly as Prior_Setup() would compute them for lmb().
## ============================================================================
set.seed(11)
n <- 25L
x1 <- rnorm(n)
x2 <- 0.7 * x1 + rnorm(n, sd = 0.6)   # correlated with x1: cor(x1,x2) ~ 0.75-0.8
X <- cbind(1, x1, x2)
beta_true <- c(1, 2, -1.5)
y <- as.numeric(X %*% beta_true + rnorm(n, sd = 1.2))

Q       <- crossprod(X)                 # base_A
beta_hat <- as.numeric(solve(Q, crossprod(X, y)))
resid   <- y - X %*% beta_hat
sigma2_hat <- sum(resid^2) / (n - 3)
V0      <- sigma2_hat * solve(Q)          # vcov(lm(y~x1+x2)) exactly
cat("cor(x1,x2) =", cor(x1, x2), "\n")
cat("V0 (vcov of full OLS fit):\n"); print(V0)

Rq <- chol(Q); Qinvhalf <- solve(Rq)

## ============================================================================
## Step 2: vector-pwt Sigma, exactly Prior_Setup()'s formula (R/prior.R:760-765)
## ============================================================================
try_pwt <- function(pwt, label, offsets_per_dim = c(-3, -2, -1, 0, 1, 2, 3)) {
  s <- sqrt((1 - pwt) / pwt)
  Sigma <- V0 * outer(s, s)
  P <- solve(Sigma); P <- 0.5 * (P + t(P))
  mu <- rep(0, 3)

  K <- t(Qinvhalf) %*% P %*% Qinvhalf
  K <- 0.5 * (K + t(K))
  eig <- eigen(K, symmetric = TRUE)
  lambda <- eig$values
  Uk <- eig$vectors
  kappa <- max(lambda) / min(lambda)

  ## Realistic dispersion bounds: mimic max_disp_perc=0.99 default off a
  ## plausible shape/rate (shape=3, rate=2*sigma2_hat as a stand-in).
  shape2 <- 5; rate2 <- 2 * sigma2_hat
  low  <- 1 / qgamma(0.995, shape2, rate2)
  upp  <- 1 / qgamma(0.005, shape2, rate2)
  t_lo <- 1 / upp; t_hi <- 1 / low

  ## Enumerate faces the way EnvelopeBuild's 3-point-per-dimension grid does:
  ## theta_bar_j = beta_hat + offset_i * omega_i for each dimension, omega_i
  ## proportional to the per-dimension prior/posterior scale. We use a modest
  ## multiplier of sqrt(diag(solve(Q))) as a stand-in for omega.
  omega <- sqrt(diag(solve(Q)))
  RSS_ML <- sum(resid^2)

  ## A(d) beta_j(d) = cbars_j - B0(d); with alpha=0, base_B0=0 here, so
  ## A(d) = P + Q/d, B0(d) = P %*% mu, beta_j(d) = solve(A(d), cbars_j - P*mu).
  rss_face_at_disp <- function(d, cbars_j) {
    A <- P + Q / d
    A <- 0.5 * (A + t(A))
    B0 <- as.numeric(P %*% mu)
    beta_j <- solve(A, cbars_j - B0)
    resid_j <- y - X %*% beta_j
    sum(resid_j^2)
  }

  grid_pts <- expand.grid(o1 = offsets_per_dim, o2 = offsets_per_dim, o3 = offsets_per_dim)
  cbars_all <- vector("list", nrow(grid_pts))
  v_all <- matrix(NA_real_, nrow(grid_pts), 3)
  rss_low_all <- numeric(nrow(grid_pts))
  for (r in seq_len(nrow(grid_pts))) {
    theta_j <- beta_hat + as.numeric(grid_pts[r, ]) * omega
    cbars_j <- as.numeric(Q %*% theta_j)              # alpha=0 here -> base_B0=0
    cbars_all[[r]] <- cbars_j
    r_star_j <- cbars_j - as.numeric(P %*% mu) - as.numeric(P %*% beta_hat)
    v_all[r, ] <- as.numeric(t(Qinvhalf) %*% r_star_j)
    rss_low_all[r] <- rss_face_at_disp(low, cbars_j)
  }

  ## RSS_min_global := min_j RSS_j(low) (RSS is decreasing in d, so its max
  ## over [low,upp] is at d=low; the global min across faces of THAT max is
  ## the real algorithm's rss_min_global). Delta_j = RSS_min_global - RSS_ML
  ## is the SAME for every face (RSS_ML/Delta do not depend on j).
  RSS_min_global <- min(rss_low_all)
  Delta <- RSS_min_global - RSS_ML
  stopifnot(Delta >= -1e-8)

  worst <- NULL
  for (r in seq_len(nrow(grid_pts))) {
    w <- as.numeric((t(Uk) %*% v_all[r, ])^2)

    ub2_lo <- ub2_red(t_lo, lambda, w, Delta)
    ub2_hi <- ub2_red(t_hi, lambda, w, Delta)
    endpoint_min <- min(ub2_lo, ub2_hi)
    ex <- find_ub2_min_by_roots(lambda, w, Delta, t_lo, t_hi)
    gap_pct <- if (abs(endpoint_min) > 1e-9) 100 * (endpoint_min - ex$ub2_min) / abs(endpoint_min) else 0

    if (is.null(worst) || gap_pct > worst$gap_pct) {
      worst <- list(row = r, w = w, Delta = Delta,
                     endpoint_min = endpoint_min, exact_min = ex$ub2_min,
                     t_star = ex$t_min, gap_pct = gap_pct)
    }
  }

  cat("\n--- ", label, " (pwt = ", paste(round(pwt, 4), collapse = ", "), ") ---\n", sep = "")
  cat("K eigenvalues:", signif(lambda, 5), " kappa(K) =", signif(kappa, 5), "\n")
  cat("RSS_ML =", signif(RSS_ML, 5), " RSS_min_global =", signif(RSS_min_global, 5),
      " Delta =", signif(Delta, 5), "\n")
  cat("Worst face found: offsets =", as.numeric(grid_pts[worst$row, ]), "\n")
  cat("  w (squared coords in K-eigenbasis):", signif(worst$w, 5), "\n")
  cat("  Delta =", signif(worst$Delta, 5), "\n")
  cat("  Endpoint-only UB2_min :", signif(worst$endpoint_min, 6), "\n")
  cat("  Exact (root) UB2_min  :", signif(worst$exact_min, 6), " at t* =", signif(worst$t_star, 5), "\n")
  cat("  Gap (endpoint overstates true min by):", round(worst$gap_pct, 3), "%\n")

  invisible(list(kappa = kappa, worst = worst))
}

cat("\n============================================================\n")
cat("Scanning vector-pwt spreads (Sigma = V0 * outer(s,s), s=sqrt((1-pwt)/pwt))\n")
cat("============================================================\n")

r1 <- try_pwt(c(0.20, 0.20, 0.20), "Equal pwt (reduces to scalar/Zellner)")
r2 <- try_pwt(c(0.30, 0.10, 0.30), "Mild spread")
r3 <- try_pwt(c(0.40, 0.05, 0.40), "Moderate spread")
r4 <- try_pwt(c(0.50, 0.02, 0.50), "Larger spread")
r5 <- try_pwt(c(0.50, 0.01, 0.20), "Larger + asymmetric spread")
r6 <- try_pwt(c(0.60, 0.005, 0.10), "Extreme spread")

cat("\n============================================================\n")
cat("Summary: kappa(K) vs worst-face gap\n")
cat("============================================================\n")
res <- list(r1, r2, r3, r4, r5, r6)
summary_df <- do.call(rbind, lapply(res, function(r) {
  data.frame(kappa = r$kappa, gap_pct = r$worst$gap_pct)
}))
print(summary_df, row.names = FALSE)
