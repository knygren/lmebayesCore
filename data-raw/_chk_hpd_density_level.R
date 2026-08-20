## For X ~ N(0, I_n):  pi(x)/pi(0) = exp(-||x||^2/2) = exp(-Xi),  Xi ~ Gamma(n/2,1).
## The 99% HPD region is {Xi <= q99}, q99 = qgamma(0.99, n/2), and the density
## threshold defining it -- the level below which only 1% of the mass sits -- is
## exp(-q99), expressed below as a percentage of the modal density.
##
## The point of the last two columns: the HPD region is a BALL of radius
## sqrt(2 q99), but any single coordinate (or any unit linear functional, by
## rotational invariance) only ranges over +-2.576 at the same 99% level.  The
## ball is far wider than any one direction needs.

tab <- do.call(rbind, lapply(c(1L, 2L, 3L, 5L, 10L, 20L, 30L, 60L, 100L), function(n) {
  q99 <- qgamma(0.99, n / 2)
  data.frame(n = n,
             q99 = q99,
             pct_of_mode = 100 * exp(-q99),
             hpd_radius = sqrt(2 * q99),
             coord_99 = qnorm(0.995),
             radius_over_coord = sqrt(2 * q99) / qnorm(0.995))
}))
cat("=== 99% HPD density threshold for N(0, I_n) ===\n")
print(tab, row.names = FALSE, digits = 4)

cat("\nRead off n = 30: the 99% HPD region reaches down to", 
    sprintf("%.2e%% of the modal density,\n", 100 * exp(-qgamma(0.99, 15))),
    "so a certified set at 1e-10 percent of the mode is not oversized -- that\n",
    "is simply what a 99% region looks like in 30 dimensions.\n")

## ---------------------------------------------------------------------------
## r(n) defined directly as asked: the level at which
##     P[ pi(beta | 0, I_n) / pi(0 | 0, I_n) <= exp(-r) ] = 0.01.
## Since the ratio is exp(-Xi) with Xi = ||beta||^2/2 ~ Gamma(n/2, 1), the event
## is {Xi >= r} and r(n) = qgamma(0.99, n/2) exactly.  Cornish-Fisher on
## Gamma(k,1) (mean k, var k, skew 2/sqrt k) gives k + z sqrt(k) + (z^2-1)/3.
zz <- qnorm(0.99)
rn <- do.call(rbind, lapply(c(1L, 2L, 3L, 4L, 5L, 8L, 10L, 15L, 20L, 30L, 40L,
                              50L, 60L, 80L, 100L, 200L, 500L, 1000L), function(n) {
  k <- n / 2
  r <- qgamma(0.99, k)
  data.frame(n = n, r_n = r, approx = k + zz * sqrt(k) + (zz^2 - 1) / 3,
             r_over_n = r / n, excess_over_half_n = r - k,
             pct_of_mode = 100 * exp(-r), log10_ratio = -r / log(10))
}))
cat("\n=== r(n): P[density ratio <= exp(-r)] = 0.01 for N(0, I_n) ===\n")
print(rn, row.names = FALSE, digits = 4)
cat("\nr(n) = qgamma(0.99, n/2).  It is essentially linear: r(n) ~ n/2 + 2.326 sqrt(n/2),\n",
    "so r/n -> 1/2 from above and the excess over n/2 grows only like sqrt(n).\n",
    "The density ratio defining the region, exp(-r(n)), therefore collapses\n",
    "geometrically in n: 3.6% at n = 1, 1e-11 at n = 30, 4e-241 at n = 1000.\n",
    "The Cornish-Fisher approximation is within 0.6% for n >= 10 and within\n",
    "0.03% for n >= 100; it is only crude at n = 1-3, where skewness dominates.\n")

## ---------------------------------------------------------------------------
## Why pct_of_mode moves even though the probability is pinned at 1%.
## The 1% is a statement about MASS; exp(-r) is a statement about HEIGHT.  The
## two decouple because the shell at radius rho has volume ~ rho^(n-1).  The
## clean way to see it: compare the density ratio at a TYPICAL draw (median of
## Xi) with the ratio at the 1% contour.  Both collapse like exp(-n/2); their
## separation grows only like sqrt(n).
cat("\n=== typical draw vs the 1% contour ===\n")
tv <- do.call(rbind, lapply(c(1L, 2L, 5L, 10L, 30L, 60L, 100L), function(n) {
  med <- qgamma(0.50, n / 2)
  r   <- qgamma(0.99, n / 2)
  data.frame(n = n, Xi_median = med, r_n = r,
             ratio_at_median = exp(-med), ratio_at_r = exp(-r),
             gap_factor = exp(r - med))
}))
print(tv, row.names = FALSE, digits = 4)
cat("\nratio_at_median is the density ratio at a TYPICAL point -- already 4e-07 at\n",
    "n = 30, with no tail behaviour involved.  The certified set only has to\n",
    "reach a further factor gap_factor = exp(r - median) beyond that, and that\n",
    "gap grows like sqrt(n), not like n.  So the low absolute level of the set\n",
    "is inherited from where the mass lives, not from the size of the budget.\n")

## ---------------------------------------------------------------------------
## Direct Monte Carlo check of the claim, with no Gamma identities used:
## draw from N(0, I_n), form the density ratio pi(x)/pi(0) exactly, and read off
## its empirical 1st percentile.  If the theory is right this matches exp(-r(n)).
set.seed(11)
N <- 2e6
cat("\n=== Monte Carlo check: empirical 1st percentile of pi(x)/pi(0) ===\n")
mc <- do.call(rbind, lapply(c(1L, 5L, 30L, 100L), function(n) {
  X   <- matrix(rnorm(N * n), N, n)
  lr  <- -0.5 * rowSums(X^2)             # log pi(x) - log pi(0), exact
  data.frame(n = n,
             mc_q01_ratio    = exp(unname(quantile(lr, 0.01))),
             theory_exp_negr = exp(-qgamma(0.99, n / 2)),
             mc_log_ratio    = unname(quantile(lr, 0.01)),
             theory_neg_r    = -qgamma(0.99, n / 2))
}))
print(mc, row.names = FALSE, digits = 5)

## The resolution of the paradox: where is the most likely radius?  The radial
## density of ||x|| is proportional to rho^(n-1) exp(-rho^2/2), peaking at
## rho = sqrt(n-1).  The density ratio THERE is already exp(-(n-1)/2).
cat("\n=== the modal radius, and the density ratio at it ===\n")
md <- do.call(rbind, lapply(c(1L, 5L, 10L, 30L, 60L, 100L), function(n) {
  rho <- sqrt(max(n - 1, 0))
  data.frame(n = n, modal_radius = rho,
             ratio_at_modal_radius = exp(-rho^2 / 2),
             ## how much probability is inside the 1%-of-mode contour
             mass_where_ratio_gt_1pct = pgamma(log(100), n / 2))
}))
print(md, row.names = FALSE, digits = 4)
cat("\nThe single MOST LIKELY radius in n = 30 already sits at 5e-07 of the modal\n",
    "density, and the region where the density exceeds 1% of its peak holds\n",
    "essentially no probability at all.  Volume beats height: the shell at\n",
    "radius rho carries volume ~ rho^(n-1), and 7.13^29 / 1^29 is about 1e+24,\n",
    "which overwhelms the 1e-11 loss in height.  So it is not that 1% of the\n",
    "MASS is unusual -- it is that ALL of the mass lives where the density is\n",
    "negligible relative to the mode.\n")

## ---------------------------------------------------------------------------
## Reparameterise the level by RADIUS rather than by log-density.  Define the
## set by  pi(beta)/pi(beta_dag) <= exp(-s^2/2),  so that s = sqrt(2 r).  For
## N(0, I_n) s is literally ||beta|| in sd units, and s(1) = 2.576 recovers the
## familiar univariate 99% z-value (whereas exp(-s^2) would give 2.576/sqrt2).
## This is the SAME nested family of sets -- a monotone change of units, not a
## new bound -- but it is the scale in which the certificate is linear.
cat("\n=== the level in radius units:  ratio <= exp(-s^2/2) ===\n")
zz99 <- qnorm(0.99); eps <- 0.01
sc <- do.call(rbind, lapply(c(1L, 2L, 5L, 10L, 30L, 60L, 100L, 500L, 1000L), function(n) {
  s_g  <- sqrt(2 * qgamma(0.99, n / 2))                      # exact Gaussian
  s_p2 <- sqrt(2 * qgamma(1 - eps / (1 + eps), n))           # Proposition 2
  data.frame(n = n, s_gauss = s_g, approx = sqrt(n) + zz99 / sqrt(2),
             s_prop2 = s_p2, p2_over_gauss = s_p2 / s_g,
             per_direction_m100 = qnorm(1 - eps / 200),
             ratio_ball_vs_dir = s_g / qnorm(1 - eps / 200))
}))
print(sc, row.names = FALSE, digits = 4)
cat("\ns(n) = sqrt(n) + 1.645 to better than 0.2% for n >= 5: the radius grows like\n",
    "sqrt(n), not like n.  The alarming collapse of exp(-r) was only the\n",
    "statement that exp(-s^2/2) with s ~ sqrt(n) is small -- a radius of 7.13 sd\n",
    "at n = 30 is not extreme.  Proposition 2 inflates the radius by a factor\n",
    "sqrt(2) = 1.414 asymptotically (1.318 at n = 30), which is far more legible\n",
    "than 'a factor 2 in r'.\n",
    "In these units the certificate is LINEAR: eta* = |eta_dag| + s / sqrt(kappa)\n",
    "and the weight floor decays like exp(-s / sqrt(kappa)).\n")

## What the ball costs relative to a per-direction budget.  In the eta scale the
## floor decays like exp(-|eta|), so the ratio of radius to per-coordinate extent
## is the exponent of the loss.
cat("\n=== cost of certifying a ball instead of each direction ===\n")
cst <- do.call(rbind, lapply(c(5L, 10L, 30L, 60L, 100L), function(n) {
  rad <- sqrt(2 * qgamma(0.99, n / 2))
  ## per-direction budget with a union bound over m linear functionals
  data.frame(n = n, ball_radius = rad,
             m10  = qnorm(1 - 0.01 / (2 * 10)),
             m100 = qnorm(1 - 0.01 / (2 * 100)),
             ratio_vs_m100 = rad / qnorm(1 - 0.01 / (2 * 100)),
             weight_factor = exp(rad - qnorm(1 - 0.01 / (2 * 100))))
}))
print(cst, row.names = FALSE, digits = 4)
cat("\nm10 / m100 are the two-sided normal quantiles after a union bound over\n",
    "10 or 100 linear functionals at total escape 0.01.  They grow like\n",
    "sqrt(2 log m) -- logarithmically -- while the ball grows like sqrt(n).\n",
    "weight_factor is exp(radius - m100): the multiplicative gain in the weight\n",
    "floor from budgeting per direction rather than over the whole ball.\n")
