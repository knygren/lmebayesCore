## Smoke: use_cpp_block1_all_chains = TRUE -> one .Call (all-chains C++ loop).

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

stopifnot(exists("two_block_block1_all_chains_cpp_export", mode = "function"))

J <- 4L
n_chains <- 2L
group_levels <- c("d", "a", "c", "b")
re_names <- c("(Intercept)", "x")
ptypes <- c("dNormal", "dIndependent_Normal_Gamma")
X_int <- matrix(1, nrow = J, ncol = 1)
X_slope <- cbind(x = c(0.1, 0.2, 0.3, 0.4), z = c(1, 0, 1, 0))
rownames(X_int) <- sort(unique(group_levels))
rownames(X_slope) <- sort(unique(group_levels))
design <- list(
  y = rbinom(J * 2L, 1, 0.3),
  Z = matrix(rnorm(J * 2L * 2L), ncol = 2L),
  groups = factor(rep(group_levels, each = 2L), levels = sort(unique(group_levels))),
  X_hyper = list("(Intercept)" = X_int, "x" = X_slope),
  re_coef_names = re_names
)
block1_prior <- list(P = diag(c(10, 20)), ddef = TRUE)
batch <- glmbayesCore:::.two_block_batch_init(
  n_chains     = n_chains,
  start_fixef  = list(
    "(Intercept)" = c("(Intercept)" = 0.5),
    "x" = c(x = 0.1, z = -0.2)
  ),
  b_start      = matrix(0, nrow = J, ncol = length(re_names),
                        dimnames = list(group_levels, re_names)),
  tau2_start   = c("(Intercept)" = 1, "x" = 2),
  re_names     = re_names,
  group_levels = group_levels
)
iters_before <- batch$iters_ranef

batch2 <- glmbayesCore:::.two_block_block1_all_chains(
  batch, design, block1_prior, stats::binomial(), ptypes,
  use_cpp_block1 = TRUE,
  use_cpp_block1_all_chains = TRUE
)
stopifnot(identical(dim(batch2$b), dim(batch$b)))
stopifnot(all(is.finite(batch2$b)))
stopifnot(length(batch2$iters_ranef) == n_chains)
stopifnot(all(is.finite(batch2$iters_ranef)))
stopifnot(any(batch2$iters_ranef >= iters_before))

message("test_block1_all_chains_cpp.R: OK")
