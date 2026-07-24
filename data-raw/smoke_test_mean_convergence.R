devtools::load_all(".", quiet = TRUE)

## --- Numeric sanity checks on the band helper itself -----------------
band_fn <- lmebayesCore:::.convergence_mean_band
n_chains_chk <- 40L
b_exact <- band_fn(n_chains_chk, 0.95, exact_ref = TRUE)
b_empir <- band_fn(n_chains_chk, 0.95, exact_ref = FALSE)
cat("exact-ref band     (normal/sqrt(n))  :", b_exact, "\n")
cat("empirical-ref band (t, wider)        :", b_empir, "\n")
stopifnot(
  isTRUE(all.equal(b_exact, c(-1, 1) * stats::qnorm(0.975) / sqrt(n_chains_chk))),
  isTRUE(all.equal(b_empir, c(-1, 1) * stats::qt(0.975, n_chains_chk - 1) / sqrt(n_chains_chk))),
  b_empir[1L] < b_exact[1L], b_empir[2L] > b_exact[2L]  ## t band strictly wider
)
cat("OK: t-band is wider than normal band, and qt(df=Inf) matches qnorm() limit\n\n")

out_dir <- tempdir()

set.seed(1)
n_chains <- 60
p1 <- 4; p2 <- 1
re_names <- c("k1", "k2")
## Chains drift geometrically toward a nonzero target so mean(l) has a
## real, shrinking bias to visualize (not just noise around 0 from sweep 1).
target <- list(k1 = c(2, -1, 0.5, 3), k2 = 4)
rate <- 0.6
mk_fixef <- function(m) {
  ## mean(m) -> target as m grows (bias target*rate^m shrinks to 0).
  list(
    k1 = matrix(
      rep(target$k1 * (1 - rate^m), each = n_chains) + rnorm(n_chains * p1, sd = 1),
      n_chains, p1, dimnames = list(NULL, paste0("x", 1:p1))
    ),
    k2 = matrix(
      rep(target$k2 * (1 - rate^m), each = n_chains) + rnorm(n_chains * p2, sd = 1),
      n_chains, p2, dimnames = list(NULL, paste0("z", 1:p2))
    )
  )
}

sweep_stats <- list()
sweep_cov <- list()
for (m in 1:10) {
  fixef <- mk_fixef(m)
  sweep_stats[[m]] <- .two_block_snapshot_fixef_stats(fixef, re_names)
  sweep_cov[[m]] <- .two_block_snapshot_fixef_cov(fixef, re_names)
}

fixef_mode <- list(
  k1 = setNames(rnorm(p1), paste0("x", 1:p1)),
  k2 = setNames(rnorm(p2), paste0("z", 1:p2))
)

hist <- .two_block_build_sweep_history(
  "main", sweep_stats, fixef_mode, re_names, sweep_cov = sweep_cov
)
cat("built hist OK; cov_by_sweep len =", length(hist$cov_by_sweep), "\n")

cat("calling plot_mean_convergence (non-whitened, empirical ref)...\n")
grDevices::png(file.path(out_dir, "mean_conv_base_nonwhitened.png"), width = 900, height = 700)
plot_mean_convergence(hist, whitened = FALSE, n_chains = n_chains)
grDevices::dev.off()
cat("non-whitened OK\n")

cat("calling plot_mean_convergence (whitened, empirical ref)...\n")
grDevices::png(file.path(out_dir, "mean_conv_base_whitened.png"), width = 900, height = 700)
plot_mean_convergence(hist, whitened = TRUE, n_chains = n_chains)
grDevices::dev.off()
cat("whitened OK\n")

## Single-series case (has_legend == FALSE path)
grDevices::png(file.path(out_dir, "mean_conv_base_singleseries.png"), width = 900, height = 700)
plot_mean_convergence(hist, whitened = FALSE, n_chains = n_chains, coef_focus = list(c("k1", "x1")))
grDevices::dev.off()
cat("single-series OK\n")

## No band (n_chains = NULL)
plot_mean_convergence(hist, whitened = FALSE)
cat("no-band OK\n")

## ggplot engine, if available
if (requireNamespace("ggplot2", quietly = TRUE)) {
  grDevices::png(file.path(out_dir, "mean_conv_ggplot_nonwhitened.png"), width = 900, height = 700)
  plot_mean_convergence(hist, whitened = FALSE, n_chains = n_chains, engine = "ggplot")
  grDevices::dev.off()
  grDevices::png(file.path(out_dir, "mean_conv_ggplot_whitened.png"), width = 900, height = 700)
  plot_mean_convergence(hist, whitened = TRUE, n_chains = n_chains, engine = "ggplot")
  grDevices::dev.off()
  cat("ggplot OK\n")
}

## --- Exact-reference path (design/measurement_prior_list) -------------
## Small synthetic known-vcov Gaussian model so lmerb_posterior_mean()/
## lmerb_posterior_covariance() both resolve.
set.seed(2)
J <- 15
n_j <- 6
grp <- factor(rep(seq_len(J), each = n_j))
D <- cbind(1, rnorm(J * n_j))
colnames(D) <- c("(Intercept)", "x1")
W_int <- cbind(1, rnorm(J))
colnames(W_int) <- c("(Intercept)", "w1")
W_slope <- matrix(1, nrow = J, ncol = 1, dimnames = list(NULL, "(Intercept)"))

b_true <- cbind(
  rnorm(J, mean = W_int %*% c(5, 2), sd = 0.3),
  rnorm(J, mean = 1, sd = 0.2)
)
y <- rowSums(D * b_true[grp, ]) + rnorm(J * n_j, sd = 0.5)

design <- list(
  y = y, Z = D, groups = grp,
  X_hyper = list("(Intercept)" = W_int, x1 = W_slope),
  re_coef_names = c("(Intercept)", "x1")
)
mpl <- list(
  dispersion_ranef = 0.5^2,
  Sigma_ranef = diag(c(0.3, 0.2))^2,
  prior_list = list(
    "(Intercept)" = list(
      mu_fixef = c(0, 0), Sigma_fixef = diag(2) * 100, dispersion_fixef = 1
    ),
    x1 = list(mu_fixef = 0, Sigma_fixef = matrix(100), dispersion_fixef = 1)
  )
)

pm <- lmerb_posterior_mean(design, mpl)
cat("exact posterior mean (fixef):\n")
print(pm$fixef)

## Build a synthetic sweep history whose chains converge toward pm$fixef,
## for the exact-reference plotting path.
re_names2 <- design$re_coef_names
target2 <- pm$fixef
mk_fixef2 <- function(m) {
  ## mean(m) -> pm$fixef (the exact posterior mean) as m grows.
  stats::setNames(lapply(re_names2, function(k) {
    p_k <- length(target2[[k]])
    matrix(
      rep(target2[[k]] * (1 - rate^m), each = n_chains) + rnorm(n_chains * p_k, sd = 1),
      n_chains, p_k, dimnames = list(NULL, names(target2[[k]]))
    )
  }), re_names2)
}
sweep_stats2 <- list(); sweep_cov2 <- list()
for (m in 1:10) {
  fixef2 <- mk_fixef2(m)
  sweep_stats2[[m]] <- .two_block_snapshot_fixef_stats(fixef2, re_names2)
  sweep_cov2[[m]] <- .two_block_snapshot_fixef_cov(fixef2, re_names2)
}
fixef_mode2 <- stats::setNames(lapply(re_names2, function(k) target2[[k]] * 0), re_names2)
hist2 <- .two_block_build_sweep_history(
  "main", sweep_stats2, fixef_mode2, re_names2, sweep_cov = sweep_cov2
)

cat("calling plot_mean_convergence (non-whitened, EXACT ref)...\n")
grDevices::png(file.path(out_dir, "mean_conv_exact_nonwhitened.png"), width = 900, height = 700)
plot_mean_convergence(
  hist2, design = design, measurement_prior_list = mpl,
  whitened = FALSE, n_chains = n_chains
)
grDevices::dev.off()
cat("exact non-whitened OK\n")

cat("calling plot_mean_convergence (whitened, EXACT ref)...\n")
grDevices::png(file.path(out_dir, "mean_conv_exact_whitened.png"), width = 900, height = 700)
plot_mean_convergence(
  hist2, design = design, measurement_prior_list = mpl,
  whitened = TRUE, n_chains = n_chains
)
grDevices::dev.off()
cat("exact whitened OK\n")

## Sanity: plot_var_convergence() still works unchanged after the rename.
cat("calling plot_var_convergence (regression check)...\n")
plot_var_convergence(hist, whitened = FALSE, n_chains = n_chains)
plot_var_convergence(hist, whitened = TRUE, n_chains = n_chains)
cat("plot_var_convergence regression OK\n")

cat("ALL OK\n")
