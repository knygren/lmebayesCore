## Trace v4 sweep-outer: n = 4 chains, verbose start/end fixef per (m, i).
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all(recompile = TRUE)

set.seed(11)
J <- 8L
n_per <- 20L
grp <- factor(rep(sprintf("g%02d", seq_len(J)), each = n_per))
w_j <- round(rnorm(J), 2)
z1 <- rnorm(J * n_per)
x_re <- cbind(`(Intercept)` = 1, slope = z1)
X_int <- cbind(1, w_j)
rownames(X_int) <- levels(grp)
x_hyper <- list(`(Intercept)` = X_int, slope = matrix(1, J, 1))
pfam_list <- list(
  `(Intercept)` = dNormal(mu = c(0, 0), Sigma = diag(4, 2), dispersion = 0.16),
  slope = dNormal(mu = 0, Sigma = diag(4, 1), dispersion = 0.16)
)
fixef_start <- list(`(Intercept)` = c(0, 0), slope = 0)
y <- rnorm(length(z1))
pl1 <- list(Sigma = diag(0.25, 2L), dispersion = 0.25)

cat("=== v4 sweep-outer trace (n = 4, m = 3) ===\n")
fit_v4 <- two_block_rNormal_reg_v4(
  n = 4L, y = y, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1, pfamily_list = pfam_list,
  fixef_start = fixef_start, m_convergence = 3L,
  progbar = FALSE, verbose = TRUE
)

cat("\n=== final fixef draws (slope col) ===\n")
print(fit_v4$fixef_draws$slope)

cat("\n=== v3 reference same settings ===\n")
fit_v3 <- two_block_rNormal_reg_v3(
  n = 4L, y = y, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1, pfamily_list = pfam_list,
  fixef_start = fixef_start, m_convergence = 3L,
  seed = 42L, seed_offset = 0L, progbar = FALSE, verbose = FALSE
)
print(fit_v3$fixef_draws$slope)
cat("max |v3 - v4| fixef:", max(abs(unlist(fit_v3$fixef_draws) - unlist(fit_v4$fixef_draws))), "\n")
