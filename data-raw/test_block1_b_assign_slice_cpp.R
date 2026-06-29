## Parity: C++ batch_b_assign_slice vs R batch$b[, , i] <- b_draw (step C).

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

stopifnot(exists(".two_block_batch_b_assign_slice_cpp", mode = "function",
                 where = asNamespace("glmbayesCore")))

J <- 4L
p_re <- 2L
n_chains <- 5L
group_levels <- c("d", "a", "c", "b")
re_names <- c("(Intercept)", "x")

set.seed(1L)
b <- array(
  rnorm(J * p_re * n_chains),
  dim = c(J, p_re, n_chains),
  dimnames = list(group_levels, re_names, NULL)
)

for (i in seq_len(n_chains)) {
  b_draw <- matrix(
    rnorm(J * p_re),
    nrow = J,
    ncol = p_re,
    dimnames = list(group_levels, re_names)
  )
  b_r <- glmbayesCore:::.two_block_batch_b_assign_slice_r(b + 0, i, b_draw)
  b_cpp <- glmbayesCore:::.two_block_batch_b_assign_slice_cpp(b + 0, i, b_draw)
  stopifnot(identical(dim(b_r), dim(b_cpp)))
  stopifnot(identical(dimnames(b_r), dimnames(b_cpp)))
  if (!identical(b_r, b_cpp)) {
    stop(
      "slice ", i, " mismatch: max abs diff = ",
      max(abs(b_r - b_cpp)),
      call. = FALSE
    )
  }
}

message("test_block1_b_assign_slice_cpp.R: OK")
