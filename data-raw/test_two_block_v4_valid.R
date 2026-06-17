if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

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
colnames(X_int) <- c("(Intercept)", "w")
X_slp <- matrix(1, J, 1L)
x_hyper <- list(`(Intercept)` = X_int, slope = X_slp)
pfam_list <- list(
  `(Intercept)` = dNormal(mu = c(0, 0), Sigma = diag(4, 2L), dispersion = 0.16),
  slope = dNormal(mu = 0, Sigma = diag(4, 1L), dispersion = 0.16)
)
fixef_start <- list(
  `(Intercept)` = c(0, 0),
  slope = 0
)
y <- rnorm(length(z1))
pl1 <- list(Sigma = diag(0.25, 2L), dispersion = 0.25)
n <- 5L
m <- 3L

fit_v4 <- two_block_rNormal_reg_v4(
  n = n, y = y, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1, pfamily_list = pfam_list,
  fixef_start = fixef_start, m_convergence = m,
  progbar = FALSE
)

stopifnot(identical(dim(fit_v4$fixef_draws[["(Intercept)"]]), c(n, 2L)))
stopifnot(all(is.finite(unlist(fit_v4$fixef_draws))))
stopifnot(all(fit_v4$iters_fixef_draws == m))
stopifnot(nrow(fit_v4$coefficients) == n * J)

cat("v4 valid: dims OK, iters = m_convergence, finite draws\n")
cat("PASS\n")
