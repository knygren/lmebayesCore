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
seed <- 42L
seed_offset <- 0L

common <- list(
  y = y, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1, pfamily_list = pfam_list,
  fixef_start = fixef_start, m_convergence = m,
  progbar = FALSE
)

fit_v3 <- do.call(two_block_rNormal_reg_v3, c(
  list(n = n, seed = seed, seed_offset = seed_offset),
  common
))

fe_ref <- matrix(NA, n, 2)
disp_ref <- matrix(NA, n, 2)
iters_ref <- matrix(NA, n, 2)
coef_ref <- vector("list", n)
for (i in seq_len(n)) {
  fit_i <- do.call(two_block_rNormal_reg_v2, c(
    list(
      n = 1L,
      seed = as.integer(seed + seed_offset + i)
    ),
    common
  ))
  fe_ref[i, ] <- fit_i$fixef_draws[["(Intercept)"]][1, ]
  disp_ref[i, ] <- fit_i$dispersion_fixef_draws[1, ]
  iters_ref[i, ] <- fit_i$iters_fixef_draws[1, ]
  coef_ref[[i]] <- fit_i$coefficients
}

fe_v3 <- fit_v3$fixef_draws[["(Intercept)"]]
max_fe <- max(abs(fe_v3 - fe_ref))
max_disp <- max(abs(fit_v3$dispersion_fixef_draws - disp_ref))
max_iters <- max(abs(fit_v3$iters_fixef_draws - iters_ref))
coef_v3 <- fit_v3$coefficients
coef_ref_rbind <- do.call(rbind, coef_ref)
re_cols <- c("(Intercept)", "slope")
max_coef <- max(abs(as.matrix(coef_v3[, re_cols, drop = FALSE]) -
                      as.matrix(coef_ref_rbind[, re_cols, drop = FALSE])))

cat("max abs diff v3 vs v2-per-chain (fixef):", max_fe, "\n")
cat("max abs diff v3 vs v2-per-chain (disp):", max_disp, "\n")
cat("max abs diff v3 vs v2-per-chain (iters):", max_iters, "\n")
cat("max abs diff v3 vs v2-per-chain (coef):", max_coef, "\n")
stopifnot(max_fe < 1e-10, max_disp < 1e-10, max_iters < 1e-10, max_coef < 1e-10)
cat("PASS\n")
