## Parity: C++ two_block_reorder_b_to_group_levels vs R reference.

pkg_root <- if (file.exists("DESCRIPTION")) "." else if (file.exists("../DESCRIPTION")) ".."
if (nzchar(pkg_root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(pkg_root, quiet = TRUE)
} else {
  library(glmbayesCore)
}

stopifnot(exists(".two_block_block1_reorder_b_r", mode = "function", where = asNamespace("glmbayesCore")))
stopifnot(exists(".two_block_reorder_b_to_group_levels_cpp", mode = "function",
                 where = asNamespace("glmbayesCore")))

group_levels <- c("c", "a", "b")
b <- matrix(
  c(1, 2, 3,
    10, 20, 30),
  nrow = 3, ncol = 2,
  byrow = TRUE,
  dimnames = list(c("a", "b", "c"), c("V1", "V2"))
)

r_out <- glmbayesCore:::.two_block_block1_reorder_b_r(b, group_levels)
cpp_out <- glmbayesCore:::.two_block_reorder_b_to_group_levels_cpp(
  b, rownames(b), group_levels
)
stopifnot(identical(r_out, cpp_out))
stopifnot(identical(rownames(r_out), group_levels))

## NULL rownames: both paths unchanged (numeric storage; C++ export is double)
b2 <- matrix(1:6, nrow = 3, ncol = 2) + 0
stopifnot(identical(
  glmbayesCore:::.two_block_block1_reorder_b_r(b2, group_levels),
  glmbayesCore:::.two_block_reorder_b_to_group_levels_cpp(b2, NULL, group_levels)
))

## Already aligned
b3 <- matrix(
  c(1, 2, 3),
  nrow = 3, ncol = 1,
  dimnames = list(group_levels, "V1")
)
stopifnot(identical(
  glmbayesCore:::.two_block_block1_reorder_b_r(b3, group_levels),
  glmbayesCore:::.two_block_reorder_b_to_group_levels_cpp(b3, rownames(b3), group_levels)
))

message("test_block1_reorder_b_cpp.R: OK")
