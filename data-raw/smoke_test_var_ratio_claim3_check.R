## Sanity check: does plot_sweep_history_var_ratio(whitened=TRUE) correctly
## recover eigenvalues <= 1 (up to small sampling noise) for DATA THAT
## ACTUALLY FOLLOWS Claim 3's growth structure (chains start coincident,
## cross-chain covariance grows monotonically toward Sigma_ref), as opposed
## to the earlier smoke test which used unstructured iid noise across
## sweeps (no real growth pattern, so ratios > 1 there are expected/normal).
devtools::load_all(".", quiet = TRUE)

set.seed(42)
n_chains <- 3000L
P <- 3L
re_names <- c("k1")
Sigma_ref_half <- diag(c(2, 1, 0.5))           ## target sd's at l = Inf
A_decay <- c(0.5, 0.8, 0.95)                    ## per-eigendirection contraction rates
L <- 10L

sweep_stats <- list(); sweep_cov <- list()
for (l in seq_len(L)) {
  sd_l <- sqrt(1 - A_decay^(2 * l))             ## Claim 3: Sigma^(l) = Sigma_ref^.5 (I-A^2l) Sigma_ref^.5
  z <- matrix(rnorm(n_chains * P), n_chains, P)
  x <- z %*% diag(sd_l) %*% Sigma_ref_half       ## n_chains x P draw with the right cross-chain covariance
  colnames(x) <- paste0("x", seq_len(P))
  fixef <- list(k1 = x)
  sweep_stats[[l]] <- .two_block_snapshot_fixef_stats(fixef, re_names)
  sweep_cov[[l]] <- .two_block_snapshot_fixef_cov(fixef, re_names)
}

fixef_mode <- list(k1 = setNames(rep(0, P), paste0("x", seq_len(P))))
hist <- .two_block_build_sweep_history("main", sweep_stats, fixef_mode, re_names, sweep_cov = sweep_cov)

grDevices::png("data-raw/claim3_check.png", width = 900, height = 700)
plot_sweep_history_var_ratio(hist, whitened = TRUE, n_chains = n_chains)
grDevices::dev.off()
cat("DONE -- if the implementation is correct, all eigenvalue traces should stay\n")
cat("at or just barely above/below 1 (small n_chains=3000 sampling noise only),\n")
cat("approaching 1 from BELOW as sweep grows, NOT persistently sitting well above 1.\n")
