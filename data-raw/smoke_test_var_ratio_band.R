devtools::load_all(".", quiet = TRUE)

## --- Numeric sanity checks on the band helper itself -----------------
band_fn <- lmebayesCore:::.sweep_var_ratio_naive_band
n_chains_chk <- 40L
b_exact <- band_fn(n_chains_chk, 0.95, exact_ref = TRUE)
b_empir <- band_fn(n_chains_chk, 0.95, exact_ref = FALSE)
cat("exact-ref band     (chi-sq/df)      :", b_exact, "\n")
cat("empirical-ref band (F, wider)       :", b_empir, "\n")
stopifnot(
  isTRUE(all.equal(b_exact, stats::qchisq(c(0.025, 0.975), n_chains_chk - 1) / (n_chains_chk - 1))),
  isTRUE(all.equal(b_empir, stats::qf(c(0.025, 0.975), n_chains_chk - 1, n_chains_chk - 1))),
  b_empir[1L] < b_exact[1L], b_empir[2L] > b_exact[2L]  ## F band strictly wider
)
cat("OK: F-band is wider than chi-sq band, and chi-sq/df matches qf(df1, Inf) limit\n\n")

## Write PNGs under the R session tempdir() (deleted on exit); change
## `out_dir` to a persistent path if you want to keep/inspect them.
out_dir <- tempdir()

set.seed(1)
n_chains <- 40
p1 <- 4; p2 <- 1
re_names <- c("k1", "k2")
mk_fixef <- function() {
  list(
    k1 = matrix(rnorm(n_chains * p1, sd = 1), n_chains, p1,
                dimnames = list(NULL, paste0("x", 1:p1))),
    k2 = matrix(rnorm(n_chains * p2, sd = 1), n_chains, p2,
                dimnames = list(NULL, paste0("z", 1:p2)))
  )
}

sweep_stats <- list()
sweep_cov <- list()
for (m in 1:8) {
  fixef <- mk_fixef()
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

grDevices::png(file.path(out_dir, "band_base_nonwhitened.png"), width = 900, height = 700)
plot_sweep_history_var_ratio(hist, whitened = FALSE, n_chains = n_chains)
grDevices::dev.off()
cat("wrote", file.path(out_dir, "band_base_nonwhitened.png"), "\n")

grDevices::png(file.path(out_dir, "band_base_whitened.png"), width = 900, height = 700)
plot_sweep_history_var_ratio(hist, whitened = TRUE, n_chains = n_chains)
grDevices::dev.off()
cat("wrote", file.path(out_dir, "band_base_whitened.png"), "\n")

## Single-series case (has_legend == FALSE path)
grDevices::png(file.path(out_dir, "band_base_singleseries.png"), width = 900, height = 700)
plot_sweep_history_var_ratio(
  hist, whitened = FALSE, n_chains = n_chains,
  coef_focus = list(c("k1", "x1"))
)
grDevices::dev.off()
cat("wrote", file.path(out_dir, "band_base_singleseries.png"), "\n")

## conf_level override
grDevices::png(file.path(out_dir, "band_base_conf99.png"), width = 900, height = 700)
plot_sweep_history_var_ratio(hist, whitened = FALSE, n_chains = n_chains, conf_level = 0.99)
grDevices::dev.off()
cat("wrote", file.path(out_dir, "band_base_conf99.png"), "\n")

## No band (n_chains = NULL) -- confirm nothing changes / no error
grDevices::png(file.path(out_dir, "band_base_noband.png"), width = 900, height = 700)
plot_sweep_history_var_ratio(hist, whitened = FALSE)
grDevices::dev.off()
cat("wrote", file.path(out_dir, "band_base_noband.png"), "\n")

## ggplot engine, if available
if (requireNamespace("ggplot2", quietly = TRUE)) {
  grDevices::png(file.path(out_dir, "band_ggplot_nonwhitened.png"), width = 900, height = 700)
  plot_sweep_history_var_ratio(hist, whitened = FALSE, n_chains = n_chains, engine = "ggplot")
  grDevices::dev.off()
  cat("wrote", file.path(out_dir, "band_ggplot_nonwhitened.png"), "\n")
}

## Error handling checks
tryCatch(
  plot_sweep_history_var_ratio(hist, n_chains = 1),
  error = function(e) cat("OK, expected error for n_chains=1:", conditionMessage(e), "\n")
)
tryCatch(
  plot_sweep_history_var_ratio(hist, n_chains = n_chains, conf_level = 1.2),
  error = function(e) cat("OK, expected error for conf_level=1.2:", conditionMessage(e), "\n")
)

cat("ALL OK\n")
