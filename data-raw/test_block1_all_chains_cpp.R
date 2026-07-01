## Smoke: .two_block_block1_all_chains_via_cpp (bulk C++ orchestrator; not production path).

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

stopifnot(exists(".two_block_block1_all_chains_via_cpp", mode = "function",
                 where = asNamespace("glmbayesCore")))

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
batch <- glmbayesCore:::.rGLMM_sweep_initialize(
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
family <- stats::binomial()

glmbayesCore:::.two_block_block1_all_chains_via_cpp(
  n            = batch$n,
  fixef        = batch$fixef,
  tau2         = batch$tau2,
  b            = batch$b,
  iters_ranef  = batch$iters_ranef,
  re_names     = batch$re_names,
  group_levels = batch$group_levels,
  design       = design,
  block1_prior = block1_prior,
  family       = family,
  ptypes       = ptypes,
  progbar      = FALSE
)

stopifnot(identical(dim(batch$b), c(J, length(re_names), n_chains)))
stopifnot(all(is.finite(batch$b)))
stopifnot(length(batch$iters_ranef) == n_chains)
stopifnot(all(is.finite(batch$iters_ranef)))
stopifnot(any(batch$iters_ranef >= iters_before))

message("test_block1_all_chains_cpp.R: OK")
