## Smoke: bulk C++ orchestrator matches R loop structure (valid outputs).

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

J <- 4L
n_chains <- 3L
group_levels <- c("d", "a", "c", "b")
re_names <- c("(Intercept)", "x")
ptypes <- c("dNormal", "dIndependent_Normal_Gamma")
X_int <- matrix(1, nrow = J, ncol = 1)
X_slope <- cbind(x = c(0.1, 0.2, 0.3, 0.4), z = c(1, 0, 1, 0))
rownames(X_int) <- sort(unique(group_levels))
rownames(X_slope) <- sort(unique(group_levels))
design <- list(
  y = c(0L, 1L, 0L, 1L, 1L, 0L, 0L, 1L),
  Z = matrix(c(0.5, -0.2, 0.1, 0.3, -0.4, 0.2, 0.6, -0.1, 0.2, 0.4, -0.3, 0.5, 0.1, -0.2, 0.3, 0.4),
             ncol = 2L, byrow = TRUE),
  groups = factor(rep(group_levels, each = 2L), levels = sort(unique(group_levels))),
  X_hyper = list("(Intercept)" = X_int, "x" = X_slope),
  re_coef_names = re_names
)
block1_prior <- list(P = diag(c(10, 20)), ddef = TRUE)
family <- stats::binomial()

init_batch <- function() {
  glmbayesCore:::.rGLMM_sweep_initialize(
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
}

run_args <- function(batch) {
  list(
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
}

batch_r <- init_batch()
batch_c <- init_batch()
out_r <- do.call(glmbayesCore:::.two_block_block1_all_chains, run_args(batch_r))
do.call(glmbayesCore:::.two_block_block1_all_chains_via_cpp, run_args(batch_c))

stopifnot(identical(dim(out_r$b), dim(batch_c$b)))
stopifnot(identical(length(out_r$iters_ranef), length(batch_c$iters_ranef)))
stopifnot(all(is.finite(out_r$b)), all(is.finite(batch_c$b)))
stopifnot(all(is.finite(out_r$iters_ranef)), all(is.finite(batch_c$iters_ranef)))

max_b_diff <- max(abs(out_r$b - batch_c$b))
stopifnot(max_b_diff < 1e-10)
stopifnot(max(abs(out_r$iters_ranef - batch_c$iters_ranef)) < 1e-10)

message(sprintf(
  "test_block1_all_chains_parity.R: OK (max |b_r - b_c| = %.2e)",
  max_b_diff
))
