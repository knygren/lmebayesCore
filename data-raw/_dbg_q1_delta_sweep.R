# q=1 synthetic limit table across delta values.
devtools::load_all(".", quiet = TRUE)

deltas <- c(0.05, 0.1, 0.25, 0.5)
kappas <- c(0.05, 0.5, 0.95)

cat("q = 1 synthetic; eps_star = (1 + kappa)^(-1/2)\n")
cat("Calibrate Pr(w Z^2 > r^2) = delta with w = kappa/(1-kappa), d = r^2/2\n\n")

for (delta in deltas) {
  cat("=== delta =", delta, "===\n")
  cat(sprintf(
    "%-7s %-9s %-7s %-7s %-9s %-9s %-9s\n",
    "kappa", "w", "r", "d", "eps_star", "eps_d", "eps"
  ))
  for (k in kappas) {
    w <- k / (1 - k)
    r1 <- .c05_r_from_delta(delta, w)
    eps_star_k <- (1 + k)^(-1 / 2)
    Q_lb <- stats::pchisq(r1$r2 / k, df = 1)
    eps_d <- exp(-r1$d) * eps_star_k
    eps <- eps_d * Q_lb
    cat(sprintf(
      "%-7.2f %-9.2f %-7.3f %-7.3f %-9.5f %-9.4e %-9.4e\n",
      k, w, r1$r, r1$d, eps_star_k, eps_d, eps
    ))
  }
  cat("\n")
}
