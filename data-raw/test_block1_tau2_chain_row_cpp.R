## Parity: C++ batch_tau2_chain_row vs R batch$tau2[i, ] (all-chains step A).

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

stopifnot(exists(".two_block_batch_tau2_chain_row_cpp", mode = "function",
                 where = asNamespace("glmbayesCore")))

re_names <- c("(Intercept)", "x")
n_chains <- 4L
tau2 <- matrix(
  c(1, 2, 3, 4,
    5, 6, 7, 8),
  nrow = n_chains,
  ncol = length(re_names),
  byrow = TRUE,
  dimnames = list(NULL, re_names)
)

for (i in seq_len(n_chains)) {
  r_row <- glmbayesCore:::.two_block_batch_tau2_chain_row_r(tau2, i)
  cpp_row <- glmbayesCore:::.two_block_batch_tau2_chain_row_cpp(tau2, i)
  stopifnot(all.equal(r_row, cpp_row, tolerance = 0))
  stopifnot(identical(names(r_row), names(cpp_row)))
}

message("test_block1_tau2_chain_row_cpp.R: OK")
