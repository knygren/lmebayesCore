## Parity: C++ batch_iters_ranef_add vs R batch$iters_ranef[i] += iters_mean (step D).

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

stopifnot(exists(".two_block_batch_iters_ranef_add_cpp", mode = "function",
                 where = asNamespace("glmbayesCore")))

n_chains <- 5L
set.seed(2L)
iters <- rnorm(n_chains)

for (i in seq_len(n_chains)) {
  delta <- runif(1L, min = 0.5, max = 20)
  it_r <- glmbayesCore:::.two_block_batch_iters_ranef_add_r(iters + 0, i, delta)
  it_cpp <- glmbayesCore:::.two_block_batch_iters_ranef_add_cpp(iters + 0, i, delta)
  if (!identical(it_r, it_cpp)) {
    stop(
      "chain ", i, " mismatch: R=", it_r[i], " C++=", it_cpp[i],
      call. = FALSE
    )
  }
}

message("test_block1_iters_ranef_add_cpp.R: OK")
