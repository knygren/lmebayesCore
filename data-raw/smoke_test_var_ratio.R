devtools::load_all(".", quiet = TRUE)
grDevices::pdf(tempfile(fileext = ".pdf"))

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
for (m in 1:5) {
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
cat("built hist OK; cov_by_sweep len =", length(hist$cov_by_sweep), "\n")
print(hist$coef_index)

cat("calling plot_sweep_history_var_ratio (non-whitened)...\n")
plot_sweep_history_var_ratio(hist, whitened = FALSE)
cat("non-whitened OK\n")

cat("calling plot_sweep_history_var_ratio (whitened)...\n")
plot_sweep_history_var_ratio(hist, whitened = TRUE)
cat("whitened OK\n")

grDevices::dev.off()
cat("ALL OK\n")
