## TEMPORARY / SCRATCH -- sanity-check .two_block_truncated_omega_moments()
## against its two degenerate limits. Not package code.

devtools::load_all(".", quiet = TRUE)

shape <- 7.3
rate  <- 4.1

## 1) Untruncated limit: omega_L -> 0, omega_U -> Inf should recover the old
##    Omega_eff / (Omega_eff^2 / shape) formula exactly.
mom_wide <- lmebayesCore:::.two_block_truncated_omega_moments(
  shape, rate, omega_L = 0, omega_U = Inf
)
Omega_eff_old   <- shape / rate
Var_eff_old     <- Omega_eff_old^2 / shape
cat(sprintf(
  "Untruncated limit: E_omega %.10f vs %.10f (old); Var_omega %.10f vs %.10f (old)\n",
  mom_wide$E_omega, Omega_eff_old, mom_wide$Var_omega, Var_eff_old
))
stopifnot(isTRUE(all.equal(mom_wide$E_omega, Omega_eff_old, tolerance = 1e-8)))
stopifnot(isTRUE(all.equal(mom_wide$Var_omega, Var_eff_old, tolerance = 1e-8)))

## 2) Point-window (degenerate truncation) limit: as omega_U -> omega_L, the
##    truncated distribution collapses to a point mass at omega_L, so
##    E_omega -> omega_L and Var_omega -> 0.
om0 <- Omega_eff_old
eps <- 1e-6
mom_point <- lmebayesCore:::.two_block_truncated_omega_moments(
  shape, rate, omega_L = om0 - eps, omega_U = om0 + eps
)
cat(sprintf(
  "Point-window limit: E_omega %.10f vs target %.10f; Var_omega %.3e (-> 0)\n",
  mom_point$E_omega, om0, mom_point$Var_omega
))
stopifnot(isTRUE(all.equal(mom_point$E_omega, om0, tolerance = 1e-3)))
stopifnot(mom_point$Var_omega < 1e-6)

## 3) A finite truncation window strictly SHRINKS the variance relative to
##    the untruncated value (per Sec 16.6's consistency check), for a window
##    that brackets the untruncated mean without being degenerate.
mom_finite <- lmebayesCore:::.two_block_truncated_omega_moments(
  shape, rate, omega_L = 0.5 * om0, omega_U = 2 * om0
)
cat(sprintf(
  "Finite window: E_omega %.6f (untrunc %.6f); Var_omega %.6f (untrunc %.6f)\n",
  mom_finite$E_omega, Omega_eff_old, mom_finite$Var_omega, Var_eff_old
))
stopifnot(mom_finite$Var_omega < Var_eff_old)

## 4) Vectorization over `rate` (draw-varying) works elementwise.
rates <- c(2, 4, 8, 16)
mom_vec <- lmebayesCore:::.two_block_truncated_omega_moments(
  shape, rates, omega_L = 0.1, omega_U = 5
)
cat("Vectorized over rate:\n")
print(data.frame(rate = rates, E_omega = mom_vec$E_omega, Var_omega = mom_vec$Var_omega, ok = mom_vec$ok))

## 5) Near-degenerate window at an extreme rate should be flagged via `ok`,
##    not silently produce garbage/NA-propagating math.
mom_deg <- lmebayesCore:::.two_block_truncated_omega_moments(
  shape, rate = 1e6, omega_L = 0.1, omega_U = 0.2
)
cat(sprintf("Near-degenerate window at large rate: dP0 = %.3e, ok = %s\n",
            mom_deg$dP0, mom_deg$ok))

cat("\nAll sanity checks passed.\n")
