## The certified level r, read as a fraction of the MODAL density.
##   Btilde_r = {beta : pitilde(beta|y) >= delta},  r = log(pitilde(bdag|y)/delta)
## so the constraint admits everything down to
##   delta / pitilde(bdag|y) = exp(-r)
## of the modal value.  In high dimension that fraction is forced to be tiny,
## and NOT because the bound is weak: it is concentration of measure.  For a
## Gaussian posterior Xi = (1/2)||delta||^2_H ~ Gamma(n/2, 1), so the TYPICAL
## draw already sits many orders of magnitude below the mode.

n <- 30L

rs <- c(3, 6, 12, 24, 25.4665, 44.216)
om <- c(0.9781, 0.9909, 0.9974, 0.9996, 0.99966, 0.999958)   # from the J=10 run
d <- data.frame(
  r            = rs,
  pct_of_mode  = 100 * exp(-rs),
  log10_ratio  = -rs / log(10),
  escape_P2    = pgamma(rs, n,     lower.tail = FALSE) / pgamma(rs, n,     lower.tail = TRUE),
  escape_gauss = pgamma(rs, n / 2, lower.tail = FALSE) / pgamma(rs, n / 2, lower.tail = TRUE),
  omega_max    = om)
cat("=== the level r as a fraction of the modal density ===\n")
print(d, row.names = FALSE, digits = 4)

cat(sprintf("\n=== where the mass actually is, n = %d ===\n", n))
cat(sprintf("Xi ~ Gamma(n/2,1): mean %.2f, median %.2f\n", n / 2, qgamma(0.5, n / 2)))
cat(sprintf("  the MEDIAN draw sits at %.3e of the modal density (%.1f orders down)\n",
            exp(-qgamma(0.5, n / 2)), qgamma(0.5, n / 2) / log(10)))
cat(sprintf("  99%% coverage needs Xi <= %.2f, i.e. %.3e of the mode\n",
            qgamma(0.99, n / 2), exp(-qgamma(0.99, n / 2))))

cat("\n=== the typical point, by dimension ===\n")
tp <- do.call(rbind, lapply(c(1L, 2L, 6L, 12L, 30L, 60L), function(nn)
  data.frame(n = nn, median_Xi = qgamma(0.5, nn / 2),
             typical_pct_of_mode = 100 * exp(-qgamma(0.5, nn / 2)),
             r_99 = qgamma(0.99, nn / 2),
             cover99_pct_of_mode = 100 * exp(-qgamma(0.99, nn / 2)))))
print(tp, row.names = FALSE, digits = 4)
cat("\nAny level-set certificate must admit the typical set, so in n = 30 it is\n",
    "FORCED below ~1e-5 percent of the modal density before the mass bound is\n",
    "even considered.  The mode is not where the posterior lives.\n")
