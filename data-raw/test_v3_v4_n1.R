pkgload::load_all()
set.seed(11)
J <- 8L
n_per <- 20L
grp <- factor(rep(sprintf("g%02d", seq_len(J)), each = n_per))
group_levels <- levels(grp)
w_j <- round(rnorm(J), 2)
z1 <- rnorm(J * n_per)
x_re <- cbind(`(Intercept)` = 1, slope = z1)
X_int <- cbind(1, w_j)
rownames(X_int) <- group_levels
x_hyper <- list(`(Intercept)` = X_int, slope = matrix(1, J, 1))
pfam_list <- list(
  `(Intercept)` = dNormal(mu = c(0, 0), Sigma = diag(4, 2), dispersion = 0.16),
  slope = dNormal(mu = 0, Sigma = diag(4, 1), dispersion = 0.16)
)
fixef_start <- list(`(Intercept)` = c(0, 0), slope = 0)
y <- rnorm(length(z1))
pl1 <- list(Sigma = diag(0.25, 2L), dispersion = 0.25)
common <- list(
  n = 1L, y = y, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1, pfamily_list = pfam_list,
  fixef_start = fixef_start, m_convergence = 3L,
  progbar = FALSE
)
f3 <- do.call(two_block_rNormal_reg_v3, c(common, list(seed = 42L, seed_offset = 0L)))
f4 <- do.call(two_block_rNormal_reg_v4, common)
max_diff <- max(abs(unlist(f3$fixef_draws) - unlist(f4$fixef_draws)))
cat("n=1 max diff:", max_diff, "\n")
cat("v3:", unlist(f3$fixef_draws), "\n")
cat("v4:", unlist(f4$fixef_draws), "\n")
