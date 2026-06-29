## Parity: C++ two_block_build_mu_all vs R build_mu_all_r.

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

stopifnot(exists("build_mu_all_r", mode = "function", where = asNamespace("glmbayesCore")))
stopifnot(exists(".two_block_build_mu_all_cpp", mode = "function",
                 where = asNamespace("glmbayesCore")))

compare_mu_all <- function(design, fixef, group_levels = NULL) {
  r_out <- build_mu_all_r(design, fixef, group_levels)
  cpp_out <- build_mu_all(design, fixef, group_levels, use_cpp = TRUE)
  stopifnot(identical(r_out$re_coef_names, cpp_out$re_coef_names))
  stopifnot(identical(r_out$group_levels, cpp_out$group_levels))
  stopifnot(identical(dimnames(r_out$mu_all), dimnames(cpp_out$mu_all)))
  stopifnot(all.equal(r_out$mu_all, cpp_out$mu_all, tolerance = 0))
  invisible(list(r = r_out, cpp = cpp_out))
}

## Synthetic: positional X_hyper rows (no rownames), permuted group_levels
J <- 4L
group_levels <- c("d", "a", "c", "b")
re_names <- c("(Intercept)", "x")
X_int <- matrix(1, nrow = J, ncol = 1)
X_slope <- cbind(
  x = c(0.1, 0.2, 0.3, 0.4),
  z = c(1, 0, 1, 0)
)
design <- list(
  X_hyper = list("(Intercept)" = X_int, "x" = X_slope),
  re_coef_names = re_names,
  groups = factor(rep(group_levels, each = 2L), levels = sort(unique(group_levels)))
)
fixef <- list(
  "(Intercept)" = c("(Intercept)" = 1.5),
  "x" = c(x = 2, z = -0.5)
)
compare_mu_all(design, fixef, group_levels)

## Synthetic: named X_hyper rows (lookup by group level)
rownames(X_int) <- sort(unique(group_levels))
rownames(X_slope) <- sort(unique(group_levels))
design$X_hyper <- list("(Intercept)" = X_int, "x" = X_slope)
compare_mu_all(design, fixef, group_levels)

## Direct export vs R reference (same inputs as build_mu_all C++ path)
x_hyper <- lapply(design$X_hyper, as.matrix)
mu_r <- build_mu_all_r(design, fixef, group_levels)$mu_all
mu_cpp <- glmbayesCore:::.two_block_build_mu_all_cpp(
  x_hyper, fixef, design$re_coef_names, group_levels
)
stopifnot(identical(dimnames(mu_r), dimnames(mu_cpp)))
stopifnot(all.equal(mu_r, mu_cpp, tolerance = 0))

message("test_block1_build_mu_all_cpp.R: OK")
